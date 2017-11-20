from __future__ import division, print_function
import matplotlib.pyplot as plt
import numpy as np
from dateutil import parser
import pytz
import jdcal
import netCDF4 as nc
import aqdlib
import datetime as dt
import xarray as xr


def load_cdf(cdf_filename, varis, squeeze_me=False):
    """Load all variables in a cdf file to a dictionary"""

    with nc.Dataset(cdf_filename, 'r') as rg:
        RAW = {}
        for var in varis:
            RAW[var] = rg[var][:]
            if squeeze_me is True:
                RAW[var] = np.squeeze(RAW[var])

        return RAW


def add_final_metadata(ds, waves=False):

    # FIXME: might not need these history lines any more
    if not waves:
        ds.attrs['history'] = 'Processed to EPIC using aqdcdf2nc.py. ' + ds.attrs['history']

        # set DELTA_T attribute for EPIC compliance.
        ds.attrs.update({'DELTA_T': ds.attrs['AQDProfileInterval']})
    else:
        ds.attrs['history'] = 'Processed to EPIC using wvscdf2nc.py. ' + ds.attrs['history']

        # set DELTA_T attribute for EPIC compliance.
        ds.attrs.update({'DELTA_T': ds.attrs['WaveInterval']})

    # print (ds['time'][0].values.astype(str))
    ds.attrs.update({'start_time': ds['time'][0].values.astype(str),
        'stop_time': ds['time'][-1].values.astype(str)})

    return ds


def add_min_max(ds):
    """
    Add minimum and maximum values to variables in NC or CDF files
    This function assumes the data are in xarray DataArrays within Datasets
    """

    exclude = list(ds.dims)
    exclude.extend(('epic_time', 'epic_time2', 'time', 'time2', 'TIM'))

    for k in ds:
        if k not in exclude:
            ds[k].attrs.update({'minimum': ds[k].min().values, 'maximum': ds[k].max().values})

    return ds


def coord_transform(vel1, vel2, vel3, heading, pitch, roll, T, cs):
    """Perform coordinate transformation to ENU"""

    N, M = np.shape(vel1)

    u = np.zeros((N,M))
    v = np.zeros((N,M))
    w = np.zeros((N,M))

    if cs == 'ENU':
        print('Data already in Earth coordinates; doing nothing')

        u = vel1
        v = vel2
        w = vel3
        # u =  vel1 * math.cos(magvar) + vel2 * math.sin(magvar);
        # v = -vel1 * math.sin(magvar) + vel2 * math.cos(magvar);
        # w = vel3;
    elif cs == 'XYZ':
        # TODO: add XYZ
        print("Data are in XYZ coordinates; transforming to Earth coordinates")
    elif cs == 'BEAM':
        print('Data are in BEAM coordinates; transforming to Earth coordinates')

        for i in range(N):
            hh = np.pi * (heading[i] - 90) / 180;
            pp = np.pi * pitch[i] / 180;
            rr = np.pi * roll[i] / 180;

            H = np.array([[ np.cos(hh), np.sin(hh), 0],
                          [-np.sin(hh), np.cos(hh), 0],
                          [ 0,          0,          1]])

            # make tilt matrix
            P = np.array([[np.cos(pp), -np.sin(pp) * np.sin(rr), -np.cos(rr) * np.sin(pp)],
                          [0,           np.cos(rr),              -np.sin(rr)],
                          [np.sin(pp),  np.sin(rr) * np.cos(pp),  np.cos(pp) * np.cos(rr)]])

            # resulting transformation matrix
            R = np.dot(np.dot(H, P), T)

            for j in range(M):
                vel = np.dot(R, np.array([vel1[i,j], vel2[i,j], vel3[i,j]]).T)
                u[i,j] = vel[0]
                v[i,j] = vel[1]
                w[i,j] = vel[2]

    return u, v, w


def set_orientation(VEL, T, metadata):
    """
    Create depth variable depending on instrument orientation
    """
    # TODO: this code seems too complicated. also should we really be modifying the trans matrix?

    N, M = np.shape(VEL['VEL1'])

    if 'Pressure_ac' in VEL:
        Wdepth = np.nanmean(VEL['Pressure_ac']) + VEL.attrs['transducer_offset_from_bottom']
    else:
        Wdepth = np.nanmean(VEL['Pressure']) + VEL.attrs['transducer_offset_from_bottom']

    blank2 = VEL.attrs['AQDBlankingDistance'] + VEL.attrs['transducer_offset_from_bottom']
    binn = VEL.attrs['bin_size']
    blank3 = VEL.attrs['transducer_offset_from_bottom'] - VEL.attrs['AQDBlankingDistance']
    binc = VEL.attrs['bin_count']

    if VEL.attrs['orientation'] == 'UP':
        print('User instructed that instrument was pointing UP')
        VEL['depth'] = xr.DataArray(np.flipud(np.linspace(Wdepth - (binn * (M - 1) + blank2 + binn), Wdepth - (blank2 + binn), num=binc)), dims=('bindist')) # need to use flipud because 1d array
    elif VEL.attrs['orientation'] == 'DOWN':
        print('User instructed that instrument was pointing DOWN')
        T[1,:] = -T[1,:]
        T[2,:] = -T[2,:]
        VEL['depth'] = xr.DataArray(np.linspace(Wdepth - blank3 + binn, Wdepth - blank3 + binn * M, num=binc),  dims=('bindist'))

    return VEL, T


def make_bin_depth(VEL, metadata):
    """Create bin_depth variable"""

    if 'Pressure_ac' in VEL:
        VEL['bin_depth'] = VEL['Pressure_ac'] - VEL['bindist']
    else:
        VEL['bin_depth'] = VEL['Pressure'] - VEL['bindist']

    return VEL


def magvar_correct(VEL, metadata):
    """Correct for magnetic declination at site"""

    if 'magnetic_variation_at_site' in metadata:
        magvardeg = metadata['magnetic_variation_at_site']
    elif 'magnetic_variation' in metadata:
        magvardeg = metadata['magnetic_variation']
    else:
        print('No magnetic variation information provided; using zero for compass correction')
        magvardeg = 0

    print('Rotating heading and horizontal velocities by %f degrees' % magvardeg)

    VEL['Heading'] = VEL['Heading'] + magvardeg
    VEL['Heading'][VEL['Heading'] >= 360] = VEL['Heading'][VEL['Heading'] >= 360] - 360
    VEL['Heading'][VEL['Heading'] < 0] = VEL['Heading'][VEL['Heading'] < 0] + 360

    vel1 = VEL['U'].copy()
    vel2 = VEL['V'].copy()

    magvar = magvardeg * np.pi / 180

    VEL['U'] =  vel1 * np.cos(magvar) + vel2 * np.sin(magvar)
    VEL['V'] = -vel1 * np.sin(magvar) + vel2 * np.cos(magvar)

    return VEL


def create_water_depth(VEL, metadata):
    """Create water_depth variable"""

    if 'initial_instrument_height' in metadata:
        if 'Pressure_ac' in VEL:
            metadata['nominal_instrument_depth'] = np.nanmean(VEL['Pressure_ac'])
            VEL['Depth'] = metadata['nominal_instrument_depth']
            wdepth = metadata['nominal_instrument_depth'] + metadata['initial_instrument_height']
            metadata['WATER_DEPTH_source'] = 'water depth = MSL from pressure sensor, atmospherically corrected'
            metadata['WATER_DEPTH_datum'] = 'MSL'
        elif 'Pressure' in VEL:
            metadata['nominal_instrument_depth'] = np.nanmean(VEL['Pressure'])
            VEL['Depth'] = metadata['nominal_instrument_depth']
            wdepth = metadata['nominal_instrument_depth'] + metadata['initial_instrument_height']
            metadata['WATER_DEPTH_source'] = 'water depth = MSL from pressure sensor'
            metadata['WATER_DEPTH_datum'] = 'MSL'
        else:
            wdepth = metadata['WATER_DEPTH']
            metadata['nominal_instrument_depth'] = metadata['WATER_DEPTH'] - metadata['initial_instrument_height']
            VEL['Depth'] = metadata['nominal_instrument_depth']
        metadata['WATER_DEPTH'] = wdepth # TODO: why is this being redefined here? Seems redundant
    elif 'nominal_instrument_depth' in metadata:
        metadata['initial_instrument_height'] = metadata['WATER_DEPTH'] - metadata['nominal_instrument_depth']
        VEL['Depth'] = metadata['nominal_instrument_depth']

    if 'initial_instrument_height' not in metadata:
        metadata['initial_instrument_height'] = 0 # TODO: do we really want to set to zero?

    return VEL, metadata


def trim_vel(VEL, metadata, waves=False):
    """Trim velocity data depending on specified method"""

    if 'Pressure_ac' in VEL:
        print('Using atmospherically corrected pressure to trim')
        WL = VEL['Pressure_ac'] + VEL.attrs['transducer_offset_from_bottom']
        P = VEL['Pressure_ac']
    else:
        # FIXME incorporate press_ ac below
        print('Using NON-atmospherically corrected pressure to trim')
        WL = VEL['Pressure'] + VEL.attrs['transducer_offset_from_bottom']
        P = VEL['Pressure']


    if 'trim_method' in metadata:
        if 'water level' in metadata['trim_method'].lower():
            if metadata['trim_method'].lower() == 'water level':
                print('Trimming using water level')
                for var in ['U', 'V', 'W', 'AGC']:
                    VEL[var] = VEL[var].where(VEL['bindist'] < P)
                VEL.attrs['history'] = 'Trimmed velocity data using water level. '+ VEL.attrs['history']
            elif metadata['trim_method'].lower() == 'water level sl':
                print('Trimming using water level and sidelobes')
                for var in ['U', 'V', 'W', 'AGC']:
                    VEL[var] = VEL[var].where(VEL['bindist'] < P * np.cos(np.deg2rad(VEL.attrs['AQDBeamAngle'])))
                VEL.attrs['history'] = 'Trimmed velocity data using water level and sidelobes. '+ VEL.attrs['history']

            # find first bin that is all bad values
            # there might be a better way to do this using xarray and named dimensions, but this works for now
            lastbin = np.argmin(np.all(np.isnan(VEL['U']), axis=0) == False)

            # this trims so there are no all-nan rows in the data
            VEL = VEL.isel(bindist=slice(0, lastbin))

            # TODO: need to add histcomment

        # TODO: add other trim methods

    return VEL


# def time_time2_to_datetime(time, time2):
#     """Create datetime array from time and time2 values"""
#
#     times = []
#
#     for t, t2 in zip(time, time2):
#         year, mon, day, frac = jdcal.jd2gcal(t - 0.5, t2/86400000)
#         hour, minute, second = day2hms(frac)
#         times.append(dt.datetime(year, mon, day, hour, minute, second, tzinfo=pytz.utc))
#
#     return np.array(times)
#
#
# def day2hms(d):
#     """Convert fractional day value into hour, minute, second"""
#
#     frachours = d * 24
#     h = np.int(np.floor(frachours))
#     fracmins = (frachours - h) * 60
#     m = np.int(np.floor(fracmins))
#     fracsecs = (fracmins - m) * 60
#     s = np.int(fracsecs)
#
#     return h, m, s
#
#
# def hms2h(h,m,s):
#     """Convert hour, minute, second to fractional hour"""
#
#     return h + m/60 + s/60/60
#
#
# def julian(t):
#     """Compute Julian date, relying heavily on jdcal package"""
#
#     y = t.year
#     m = t.month
#     d = t.day
#     h = hms2h(t.hour, t.minute, t.second)
#     return sum(jdcal.gcal2jd(y,m,d)) + h/24 + 0.5
