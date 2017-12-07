from __future__ import division, print_function
import matplotlib.pyplot as plt
import numpy as np
from dateutil import parser
import pytz
import jdcal
import netCDF4 as nc
import datetime as dt
import xarray as xr


def add_final_aqd_metadata(ds, waves=False):

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


def set_orientation(VEL, T):
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


def make_bin_depth(VEL):
    """Create bin_depth variable"""

    if 'Pressure_ac' in VEL:
        VEL['bin_depth'] = VEL['Pressure_ac'] - VEL['bindist']
    else:
        VEL['bin_depth'] = VEL['Pressure'] - VEL['bindist']

    return VEL


def magvar_correct(ds):
    """Correct for magnetic declination at site"""

    if 'magnetic_variation_at_site' in ds.attrs:
        magvardeg = ds.attrs['magnetic_variation_at_site']
    elif 'magnetic_variation' in ds.attrs:
        magvardeg = ds.attrs['magnetic_variation']
    else:
        print('No magnetic variation information provided; using zero for compass correction')
        magvardeg = 0

    print('Rotating heading and horizontal velocities by %f degrees' % magvardeg)

    ds['Heading'] = ds['Heading'] + magvardeg
    ds['Heading'][ds['Heading'] >= 360] = ds['Heading'][ds['Heading'] >= 360] - 360
    ds['Heading'][ds['Heading'] < 0] = ds['Heading'][ds['Heading'] < 0] + 360

    vel1 = ds['U'].copy()
    vel2 = ds['V'].copy()

    magvar = magvardeg * np.pi / 180

    ds['U'] =  vel1 * np.cos(magvar) + vel2 * np.sin(magvar)
    ds['V'] = -vel1 * np.sin(magvar) + vel2 * np.cos(magvar)

    return ds


def trim_vel(ds, waves=False):
    """Trim velocity data depending on specified method"""

    if 'Pressure_ac' in ds:
        print('Using atmospherically corrected pressure to trim')
        WL = ds['Pressure_ac'] + ds.attrs['transducer_offset_from_bottom']
        P = ds['Pressure_ac']
    else:
        # FIXME incorporate press_ ac below
        print('Using NON-atmospherically corrected pressure to trim')
        WL = ds['Pressure'] + ds.attrs['transducer_offset_from_bottom']
        P = ds['Pressure']


    if 'trim_method' in ds.attrs:
        if 'water level' in ds.attrs['trim_method'].lower():
            if ds.attrs['trim_method'].lower() == 'water level':
                print('Trimming using water level')
                for var in ['U', 'V', 'W', 'AGC']:
                    ds[var] = ds[var].where(ds['bindist'] < P)
                ds.attrs['history'] = 'Trimmed velocity data using water level. '+ ds.attrs['history']
            elif ds.attrs['trim_method'].lower() == 'water level sl':
                print('Trimming using water level and sidelobes')
                for var in ['U', 'V', 'W', 'AGC']:
                    ds[var] = ds[var].where(ds['bindist'] < P * np.cos(np.deg2rad(ds.attrs['AQDBeamAngle'])))
                ds.attrs['history'] = 'Trimmed velocity data using water level and sidelobes. '+ ds.attrs['history']

            # find first bin that is all bad values
            # there might be a better way to do this using xarray and named dimensions, but this works for now
            lastbin = np.argmin(np.all(np.isnan(ds['U']), axis=0) == False)

            # this trims so there are no all-nan rows in the data
            ds = ds.isel(bindist=slice(0, lastbin))

            # TODO: need to add histcomment

        # TODO: add other trim methods

    return ds
