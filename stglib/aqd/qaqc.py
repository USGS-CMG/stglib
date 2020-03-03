from __future__ import division, print_function
import numpy as np
import math
import xarray as xr
from ..core import utils


def ds_rename(ds, waves=False):
    """
    Rename DataArrays within Dataset for EPIC compliance
    """

    varnames = {'Pressure': 'P_1',
                'Temperature': 'Tx_1211',
                'Heading': 'Hdg_1215',
                'Pitch': 'Ptch_1216',
                'Roll': 'Roll_1217',
                'Battery': 'Bat_106'}

    if 'Pressure_ac' in ds:
        varnames['Pressure_ac'] = 'P_1ac'

    if not waves:
        varnames.update(
            {'U': 'u_1205',
             'V': 'v_1206',
             'W': 'w_1204',
             'AGC': 'AGC_1202'})
    elif waves:
        varnames.update(
            {'VEL1': 'vel1_1277',
             'VEL2': 'vel2_1278',
             'VEL3': 'vel3_1279',
             'AMP1': 'AGC1_1221',
             'AMP2': 'AGC2_1222',
             'AMP3': 'AGC3_1223'})

    ds = ds.rename(varnames)

    if waves:
        for v in ['vel1_1277',
                  'vel2_1278',
                  'vel3_1279',
                  'AGC1_1221',
                  'AGC2_1222',
                  'AGC3_1223']:
            ds[v] = ds[v].expand_dims('depth', axis=-1)

    for v in ['avgamp1',
              'avgamp2',
              'avgamp3',
              'U',
              'V',
              'W',
              'Depth',
              'water_depth',
              'cellpos']:
        if v in ds:
            ds = ds.drop(v)

    return ds


def load_cdf(cdf_filename, atmpres=False):
    """
    Load raw .cdf file and, optionally, an atmospheric pressure .cdf file
    """

    ds = xr.open_dataset(cdf_filename).load()
    ds.close()

    if atmpres is not False:
        p = xr.open_dataset(atmpres).load()
        p.close()
        # TODO: check to make sure this data looks OK
        ds['Pressure_ac'] = xr.DataArray(ds['Pressure'] -
                                         p['atmpres'] -
                                         p['atmpres'].offset)

    return ds


def date_parser(year, month, day, hour, minute, second):
    """read csv and parse dates
    https://stackoverflow.com/q/27112591
    """

    return (year + '-' + month + '-' + day + ' ' +
            hour + ':' + minute + ':' + second)


def add_delta_t(ds, waves=False):
    """
    set DELTA_T attribute for EPIC compliance.
    """
    if not waves:
        ds.attrs.update({'DELTA_T': ds.attrs['AQDProfileInterval']})
    else:
        ds.attrs.update({'DELTA_T': ds.attrs['WaveInterval']})

    return ds


def make_tilt(p, r):
    return np.array([
        [math.cos(p), -math.sin(p)*math.sin(r), -math.cos(r)*math.sin(p)],
        [0,           math.cos(r),              -math.sin(r)],
        [math.sin(p),  math.sin(r)*math.cos(p),  math.cos(p)*math.cos(r)]])


def coord_transform(vel1, vel2, vel3, heading, pitch, roll, T, T_orig, cs):
    """Perform coordinate transformation to ENU"""

    N, M = np.shape(vel1)

    u = np.zeros((N, M))
    v = np.zeros((N, M))
    w = np.zeros((N, M))

    if cs == 'ENU':
        print('Data already in Earth coordinates; doing nothing')

        u = vel1
        v = vel2
        w = vel3

    elif cs == 'XYZ' or cs == 'BEAM':
        print('Data are in %s coordinates; transforming to Earth '
              'coordinates' % cs)

        for i in range(N):
            hh = np.pi * (heading[i] - 90) / 180
            pp = np.pi * pitch[i] / 180
            rr = np.pi * roll[i] / 180

            H = np.array([[math.cos(hh),  math.sin(hh), 0],
                          [-math.sin(hh), math.cos(hh), 0],
                          [0,             0,            1]])

            # make tilt matrix
            P = make_tilt(pp, rr)

            # resulting transformation matrix
            R = np.dot(np.dot(H, P), T)
            if cs == 'XYZ':
                for j in range(M):
                    vel = np.dot(
                        np.dot(R, np.linalg.inv(T_orig)),
                        np.array([vel1[i, j], vel2[i, j], vel3[i, j]]).T)
                    u[i, j] = vel[0]
                    v[i, j] = vel[1]
                    w[i, j] = vel[2]

            elif cs == 'BEAM':
                for j in range(M):
                    vel = np.dot(
                        R, np.array([vel1[i, j], vel2[i, j], vel3[i, j]]).T)
                    u[i, j] = vel[0]
                    v[i, j] = vel[1]
                    w[i, j] = vel[2]

    return u, v, w


def swap_bindist_to_depth(ds):
    return ds.swap_dims({'bindist': 'depth'})


def set_orientation(VEL, T):
    """
    Create depth variable depending on instrument orientation
    """

    if 'Pressure_ac' in VEL:
        presvar = 'Pressure_ac'
    else:
        presvar = 'Pressure'

    if 'NAVD88_ref' in VEL.attrs:
        Wdepth = (-VEL.attrs['NAVD88_ref'] -
                  VEL.attrs['transducer_offset_from_bottom'])
    else:
        Wdepth = np.nanmean(VEL[presvar])

    T_orig = T.copy()

    if VEL.attrs['orientation'] == 'UP':
        print('User instructed that instrument was pointing UP')
        # try a simpler approach
        VEL['depth'] = xr.DataArray(Wdepth - VEL['bindist'], dims='bindist')
    elif VEL.attrs['orientation'] == 'DOWN':
        print('User instructed that instrument was pointing DOWN')
        T[1, :] = -T[1, :]
        T[2, :] = -T[2, :]
        # try a simpler approach
        VEL['depth'] = xr.DataArray(Wdepth + VEL['bindist'], dims='bindist')

    return VEL, T, T_orig


def make_bin_depth(VEL, waves=False):
    """Create bin_depth variable"""

    if 'Pressure_ac' in VEL:
        pres = 'Pressure_ac'
    else:
        pres = 'Pressure'

    if not waves:
        VEL['bin_depth'] = VEL[pres] - VEL['bindist']
    else:
        VEL['bin_depth'] = VEL[pres].mean(dim='sample') - VEL['bindist']

    return VEL


def magvar_correct(ds):
    """Correct for magnetic declination at site"""

    if 'magnetic_variation_at_site' in ds.attrs:
        magvardeg = ds.attrs['magnetic_variation_at_site']
    elif 'magnetic_variation' in ds.attrs:
        magvardeg = ds.attrs['magnetic_variation']
    else:
        print('No magnetic variation information provided; '
              'using zero for compass correction')
        magvardeg = 0

    print('Rotating heading and horizontal velocities by %f degrees' %
          magvardeg)

    ds['Heading'] = ds['Heading'] + magvardeg
    ds['Heading'][ds['Heading'] >= 360] = (
        ds['Heading'][ds['Heading'] >= 360] - 360)
    ds['Heading'][ds['Heading'] < 0] = ds['Heading'][ds['Heading'] < 0] + 360

    ds['U'], ds['V'] = rotate(ds['U'],  ds['V'], magvardeg)

    return ds


def rotate(u, v, deg):
    rad = np.deg2rad(deg)
    urot = u * np.cos(rad) + v * np.sin(rad)
    vrot = -u * np.sin(rad) + v * np.cos(rad)

    return urot, vrot


def trim_vel(ds, waves=False, data_vars=['U', 'V', 'W', 'AGC']):
    """Trim velocity data depending on specified method

    Parameters
    ----------
    ds : xarray.Dataset
        The xarray Dataset
    waves: bool, optional
        Flag to determine whether these are waves data. Default False.
    data_vars : array_like
        List of variables to trim. Default ['U', 'V', 'W', 'AGC'].

    Returns
    -------
    xarray.Dataset
        Dataset with trimmed data
    """

    if ('trim_method' in ds.attrs and
            ds.attrs['trim_method'].lower() != 'none' and
            ds.attrs['trim_method'] is not None):

        if 'Pressure_ac' in ds:
            print('Using atmospherically corrected pressure to trim')
            WL = ds['Pressure_ac'] + ds.attrs['transducer_offset_from_bottom']
            P = ds['Pressure_ac']
        elif 'Pressure' in ds:
            # FIXME incorporate press_ ac below
            print('Using NON-atmospherically corrected pressure to trim')
            WL = ds['Pressure'] + ds.attrs['transducer_offset_from_bottom']
            P = ds['Pressure']

        if ds.attrs['trim_method'].lower() == 'water level':
            print('Trimming using water level')
            for var in data_vars:
                ds[var] = ds[var].where(ds['bindist'] < P)
            ds.attrs['history'] = (
                'Trimmed velocity data using water level. ' +
                ds.attrs['history'])
        elif ds.attrs['trim_method'].lower() == 'water level sl':
            print('Trimming using water level and sidelobes')
            for var in data_vars:
                ds[var] = ds[var].where(
                    ds['bindist'] <
                    P * np.cos(np.deg2rad(ds.attrs['AQDBeamAngle'])))
            ds.attrs['history'] = (
                'Trimmed velocity data using water level and sidelobes. ' +
                ds.attrs['history'])
        elif ds.attrs['trim_method'].lower() == 'bin range':
            print('Trimming using good_bins of %s' %
                  str(ds.attrs['good_bins']))
            for var in data_vars:
                ds[var] = ds[var].isel(bindist=slice(ds.attrs['good_bins'][0],
                                                     ds.attrs['good_bins'][1]))

        # find first bin that is all bad values
        # there might be a better way to do this using xarray and named
        # dimensions, but this works for now
        lastbin = np.argmin(
            np.all(np.isnan(ds[data_vars[0]].values), axis=0) == False)

        # this trims so there are no all-nan rows in the data
        ds = ds.isel(bindist=slice(0, lastbin))

        # TODO: need to add histcomment

        # TODO: add other trim methods
    else:
        print('Did not trim velocity data')

    return ds


def read_aqd_hdr(basefile):
    """
    Get instrument metadata from .hdr file
    Was formerly readAQDprfHeader.m
    """

    hdrFile = basefile + '.hdr'

    # read in whole file first; shouldn't be too bad since it's small, to
    # check and make sure it's not an HR file.
    # FIXME: this will need to be removed when we start supporting HR
    # check for HR by seeing if extended velocity range is in .hdr
    with open(hdrFile) as f:
        if 'Extended velocity range' in f.read():
            raise NotImplementedError(
                'stglib does not currently support Aquadopp HR datasets')

    f = open(hdrFile, 'r')
    row = ''

    Instmeta = {}

    while 'Hardware configuration' not in row:
        row = f.readline().rstrip()
        if 'Profile interval' in row:
            idx = row.find(' sec')
            Instmeta['AQDProfileInterval'] = int(row[38:idx])
        elif 'Number of cells' in row:
            Instmeta['AQDNumberOfCells'] = int(row[38:])
        # required here to differentiate from the wave cell size
        elif row.find('Cell size', 0, 9) != -1:
            idx = row.find(' cm')
            Instmeta['AQDCellSize'] = int(row[38:idx])
        elif 'Average interval' in row:
            idx = row.find(' sec')
            Instmeta['AQDAverageInterval'] = int(row[38:idx])
        elif 'Measurement load' in row:
            idx = row.find(' %')
            Instmeta['AQDMeasurementLoad'] = int(row[38:idx])
        elif 'Transmit pulse length' in row:
            idx = row.find(' m')
            Instmeta['AQDTransmitPulseLength'] = float(row[38:idx])
        elif 'Blanking distance' in row:
            idx = row.find(' m')
            Instmeta['AQDBlankingDistance'] = float(row[38:idx])
        elif 'Compass update rate' in row:
            idx = row.find(' sec')
            Instmeta['AQDCompassUpdateRate'] = int(row[38:idx])
        elif 'Wave measurements' in row:
            Instmeta['WaveMeasurements'] = row[38:]
        elif 'Wave - Powerlevel' in row:
            Instmeta['WavePower'] = row[38:]
        elif 'Wave - Interval' in row:
            idx = row.find(' sec')
            Instmeta['WaveInterval'] = int(row[38:idx])
        elif 'Wave - Number of samples' in row:
            Instmeta['WaveNumberOfSamples'] = int(row[38:])
        elif 'Wave - Sampling rate' in row:
            Instmeta['WaveSampleRate'] = row[38:]
        elif 'Wave - Cell size' in row:
            idx = row.find(' m')
            Instmeta['WaveCellSize'] = float(row[38:idx])
        elif 'Analog input 1' in row:
            Instmeta['AQDAnalogInput1'] = row[38:]
        elif 'Analog input 2' in row:
            Instmeta['AQDAnalogInput2'] = row[38:]
        elif 'Power output' in row:
            Instmeta['AQDAnalogPowerOutput'] = row[38:]
        elif 'Powerlevel' in row:
            Instmeta['AQDPowerLevel'] = row[38:]
        elif 'Coordinate system' in row:
            Instmeta['AQDCoordinateSystem'] = row[38:]
        elif 'Sound speed' in row:
            Instmeta['AQDSoundSpeed'] = row[38:]
        elif 'Salinity' in row:
            Instmeta['AQDSalinity'] = row[38:]
        elif 'Number of beams' in row:
            Instmeta['AQDNumberOfBeams'] = int(row[38:])
        elif 'Number of pings per burst' in row:
            Instmeta['AQDNumberOfPingsPerBurst'] = int(row[38:])
        elif 'Software version' in row:
            Instmeta['AQDSoftwareVersion'] = row[38:]
        elif 'Deployment name' in row:
            Instmeta['AQDDeploymentName'] = row[38:]
        elif 'Deployment time' in row:
            Instmeta['AQDDeploymentTime'] = row[38:]
        elif 'Comments' in row:
            Instmeta['AQDComments'] = row[38:]

    while 'Head configuration' not in row:
        row = f.readline().rstrip()
        if 'Serial number' in row:
            Instmeta['AQDSerial_Number'] = row[38:]
        elif 'Hardware revision' in row:
            Instmeta['AQDHardwareRevision'] = row[38:]
        elif 'Revision number' in row:
            Instmeta['AQDRevisionNumber'] = row[38:]
        elif 'Recorder size' in row:
            Instmeta['AQDRecorderSize'] = row[38:]
        elif 'Firmware version' in row:
            Instmeta['AQDFirmwareVersion'] = row[38:]
        elif 'Velocity range' in row:
            Instmeta['AQDVelocityRange'] = row[38:]
        elif 'Power output' in row:
            Instmeta['AQDAnalogPowerOutput'] = row[38:]
        elif 'Analog input #1 calibration (a0, a1)' in row:
            Instmeta['AQDAnalogInputCal1'] = row[38:]
        elif 'Analog input #2 calibration (a0, a1)' in row:
            Instmeta['AQDAnalogInputCal2'] = row[38:]
        elif 'Sync signal data out delay' in row:
            Instmeta['AQDSyncOutDelay'] = row[38:]
        elif 'Sync signal power down delay' in row:
            Instmeta['AQDSyncPowerDelay'] = row[38:]

    while 'Current profile cell center distance from head (m)' not in row:
        row = f.readline().rstrip()
        if 'Pressure sensor' in row:
            Instmeta['AQDPressureSensor'] = row[38:]
        elif 'Compass' in row:
            Instmeta['AQDCompass'] = row[38:]
        elif 'Tilt sensor' in row:
            Instmeta['AQDTilt'] = row[38:]
        elif 'Head frequency' in row:
            idx = row.find(' kHz')
            Instmeta['AQDFrequency'] = int(row[38:idx])
        elif 'Number of beams' in row:
            Instmeta['AQDNumBeams'] = int(row[38:])
        elif 'Serial number' in row:
            Instmeta['AQDHeadSerialNumber'] = row[38:]
        elif 'Transformation matrix' in row:
            Instmeta['AQDTransMatrix'] = np.zeros((3, 3))
            for n in np.arange(3):
                Instmeta['AQDTransMatrix'][n, :] = (
                    [float(x) for x in row[38:].split()])
                row = f.readline().rstrip()
        elif 'Pressure sensor calibration' in row:
            Instmeta['AQDPressureCal'] = row[38:]

    bd = []
    while 'Data file format' not in row:
        row = f.readline().rstrip()
        # avoid the header rule line
        if '-' not in row and row != '' and row != 'Data file format':
            bd.append(float(row.split()[1]))

    Instmeta['AQDCCD'] = np.array(bd)  # CCD = Cell Center Distance

    # infer some things based on the Aquadopp brochure
    if Instmeta['AQDFrequency'] == 400:
        Instmeta['AQDBeamWidth'] = 3.7
    elif Instmeta['AQDFrequency'] == 600:
        Instmeta['AQDBeamWidth'] = 3.0
    elif Instmeta['AQDFrequency'] == 1000:
        Instmeta['AQDBeamWidth'] = 3.4
    elif Instmeta['AQDFrequency'] == 2000:
        Instmeta['AQDBeamWidth'] = 1.7
    else:
        Instmeta['AQDBeamWidth'] = np.nan

    Instmeta['AQDBeamPattern'] = 'convex'
    Instmeta['AQDBeamAngle'] = 25

    f.close()

    return Instmeta


def check_attrs(ds, waves=False):

    # Add some metadata originally in the run scripts
    ds.attrs['nominal_sensor_depth_note'] = ('WATER_DEPTH - '
                                             'initial_instrument_height')
    ds.attrs['nominal_sensor_depth'] = (
        ds.attrs['WATER_DEPTH'] -
        ds.attrs['initial_instrument_height'])
    ds.attrs['transducer_offset_from_bottom'] = (
        ds.attrs['initial_instrument_height'])

    if ('initial_instrument_height' not in ds.attrs or
            np.isnan(ds.attrs['initial_instrument_height'])):
        ds.attrs['initial_instrument_height'] = 0

    ds.attrs['serial_number'] = ds.attrs['AQDSerial_Number']

    # update metadata from Aquadopp header to CMG standard so that various
    # profilers have the same attribute wording.  Redundant, but necessary
    if not waves:
        ds.attrs['bin_count'] = ds.attrs['AQDNumberOfCells']
        ds.attrs['bin_size'] = ds.attrs['AQDCellSize'] / 100  # from cm to m
        ds.attrs['blanking_distance'] = ds.attrs['AQDBlankingDistance']
        # Nortek lists the distance to the center of the first bin as the
        # blanking distance plus one cell size
        ds.attrs['center_first_bin'] = (ds.attrs['blanking_distance'] +
                                        ds.attrs['bin_size'])  # in m
    else:
        ds.attrs['bin_count'] = 1  # only 1 wave bin
        ds.attrs['bin_size'] = ds.attrs['WaveCellSize']
        ds.attrs['blanking_distance'] = ds.attrs['AQDBlankingDistance']
        # TODO: need to set center_first_bin after return in main
        # calling function

    ds.attrs['salinity_set_by_user'] = ds.attrs['AQDSalinity']
    ds.attrs['salinity_set_by_user_units'] = 'ppt'

    ds.attrs['frequency'] = ds.attrs['AQDFrequency']
    ds.attrs['beam_width'] = ds.attrs['AQDBeamWidth']
    ds.attrs['beam_pattern'] = ds.attrs['AQDBeamPattern']
    ds.attrs['beam_angle'] = ds.attrs['AQDBeamAngle']

    # TODO: Ideally we want to read the error, status, and orientation from
    # the .SEN file, but this requires reading the first good burst, which
    # is unknown at this point.

    ds.attrs['INST_TYPE'] = 'Nortek Aquadopp Profiler'

    return ds


def check_orientation(ds, waves=False):
    """Check instrument orientation and create variables that depend on this"""

    print('Insrument orientation:', ds.attrs['orientation'])
    print('Center_first_bin = %f' % ds.attrs['center_first_bin'])
    print('bin_size = %f' % ds.attrs['bin_size'])
    print('bin_count = %f' % ds.attrs['bin_count'])
    # TODO: these values are already in the HDR file...
    if not waves:
        bindist = np.linspace(ds.attrs['center_first_bin'],
                              (ds.attrs['center_first_bin'] +
                               ((ds.attrs['bin_count'] - 1) *
                               ds.attrs['bin_size'])),
                              num=ds.attrs['bin_count'])
    else:
        bindist = ds['cellpos'][0]

    if ds.attrs['orientation'] == 'UP':
        print('User instructed that instrument was pointing UP')
        # depth, or distance below surface, is a positive number below the
        # surface, negative above the surface, for CMG purposes and consistency
        # with ADCP
        depth = (ds.attrs['WATER_DEPTH'] -
                 ds.attrs['transducer_offset_from_bottom']) - bindist
        # TODO: this is never used
        Depth_NOTE = ('user reports uplooking bin depths = water_depth - '
                      'transducer offset from bottom - bindist')
    elif ds.attrs['orientation'] == 'DOWN':
        print('User instructed that instrument was pointing DOWN')
        depth = (ds.attrs['WATER_DEPTH'] -
                 ds.attrs['transducer_offset_from_bottom']) + bindist
        # TODO: this is never used
        Depth_NOTE = ('user reports downlooking bin depths = water_depth - '
                      'transducer_offset_from_bottom + bindist')

    if not waves:
        ds['bindist'] = xr.DataArray(bindist, dims=('bindist'), name='bindist')
        ds['depth'] = xr.DataArray(depth, dims=('bindist'), name='depth')
    else:
        ds['bindist'] = xr.DataArray([bindist],
                                     dims=('bindist'),
                                     name='bindist')
        ds['depth'] = xr.DataArray([depth], dims=('bindist'), name='depth')

    return ds


def update_attrs(ds, waves=False):
    """Define dimensions and variables in NetCDF file"""

    ds['lat'] = xr.DataArray([ds.attrs['latitude']], dims=('lat'), name='lat')
    ds['lon'] = xr.DataArray([ds.attrs['longitude']], dims=('lon'), name='lon')

    ds['TransMatrix'] = xr.DataArray(ds.attrs['AQDTransMatrix'],
                                     dims=('Tmatrix', 'Tmatrix'),
                                     name='TransMatrix')
    # Need to remove AQDTransMatrix from attrs for netCDF3 compliance
    ds.attrs.pop('AQDTransMatrix')

    ds['time'].attrs.update(
        {'standard_name': 'time',
         'axis': 'T'})

    ds['lat'].attrs.update(
        {'units': 'degree_north',
         'long_name': 'Latitude',
         'standard_name': 'latitude',
         'epic_code': 500})

    ds['lon'].attrs.update(
        {'units': 'degree_east',
         'long_name': 'Longitude',
         'standard_name': 'longitude',
         'epic_code': 502})

    if 'position_datum' in ds.attrs:
        ds['lat'].attrs['datum'] = ds.attrs['position_datum']
        ds['lon'].attrs['datum'] = ds.attrs['position_datum']

    ds['bindist'].attrs.update(
        {'units': 'm',
         'long_name': 'distance from transducer head',
         'bin_size': ds.attrs['bin_size'],
         'center_first_bin': ds.attrs['center_first_bin'],
         'bin_count': ds.attrs['bin_count'],
         'transducer_offset_from_bottom':
            ds.attrs['transducer_offset_from_bottom']})

    ds['Temperature'].attrs.update(
        {'units': 'C',
         'long_name': 'Temperature',
         'generic_name': 'temp'})

    ds['Pressure'].attrs.update(
        {'units': 'dbar',
         'long_name': 'Pressure',
         'generic_name': 'press',
         'note': ('Raw pressure from instrument, not corrected for changes '
                  'in atmospheric pressure')})

    for n in [1, 2, 3]:
        ds['VEL' + str(n)].attrs.update({
            'units': 'cm/s',
            'Type': 'scalar',
            'transducer_offset_from_bottom':
                ds.attrs['transducer_offset_from_bottom']})
        ds['AMP' + str(n)].attrs.update({
            'long_name': 'Beam ' + str(n) + ' Echo Amplitude',
            'units': 'counts',
            'Type': 'scalar',
            'transducer_offset_from_bottom':
                ds.attrs['transducer_offset_from_bottom']})

    if not waves:
        veltxt = 'current velocity'
    else:
        veltxt = 'wave-burst velocity'

    if ds.attrs['AQDCoordinateSystem'] == 'ENU':
        ds['VEL1'].attrs['long_name'] = 'Eastward ' + veltxt
        ds['VEL2'].attrs['long_name'] = 'Northward ' + veltxt
        ds['VEL3'].attrs['long_name'] = 'Vertical ' + veltxt
    elif ds.attrs['AQDCoordinateSystem'] == 'XYZ':
        ds['VEL1'].attrs['long_name'] = veltxt.capitalize() + ' in X Direction'
        ds['VEL2'].attrs['long_name'] = veltxt.capitalize() + ' in Y Direction'
        ds['VEL3'].attrs['long_name'] = veltxt.capitalize() + ' in Z Direction'
    elif ds.attrs['AQDCoordinateSystem'] == 'BEAM':
        ds['VEL1'].attrs['long_name'] = 'Beam 1 ' + veltxt
        ds['VEL2'].attrs['long_name'] = 'Beam 2 ' + veltxt
        ds['VEL3'].attrs['long_name'] = 'Beam 3 ' + veltxt

    ds['Battery'].attrs.update(
        {'units': 'Volts',
         'long_name': 'Battery Voltage'})

    ds['Pitch'].attrs.update(
        {'units': 'degrees',
         'long_name': 'Instrument Pitch'})

    ds['Roll'].attrs.update(
        {'units': 'degrees',
         'long_name': 'Instrument Roll'})

    ds['Heading'].attrs.update(
        {'units': 'degrees',
         'long_name': 'Instrument Heading',
         'datum': 'magnetic north'})

    ds['depth'].attrs.update(
        {'units': 'm',
         'long_name': 'mean water depth',
         'bin_size': ds.attrs['bin_size'],
         'center_first_bin': ds.attrs['center_first_bin'],
         'bin_count': ds.attrs['bin_count'],
         'transducer_offset_from_bottom':
            ds.attrs['transducer_offset_from_bottom']})

    ds['TransMatrix'].attrs['long_name'] = ('Transformation Matrix '
                                            'for this Aquadopp')
    if waves:
        ds['burst'].attrs.update(
            {'units': 'count',
             'long_name': 'Record number'})
        # ds['burst'].encoding['_FillValue'] = 1e35 # don't want this to have a _FillValue

        # ds['cellpos'].encoding['_FillValue'] = 1e35 # don't want this to have a _FillValue

    return ds


def ds_add_attrs(ds, waves=False):
    """
    add EPIC and CMG attributes to xarray Dataset
    """

    def add_vel_attributes(vel, dsattrs):
        vel.attrs.update(
            {'units': 'cm/s',
             'data_cmnt': ('Velocity in shallowest bin is often suspect and '
                           'should be used with caution')})

        # TODO: why do we only do trim_method for Water Level SL?
        if ('trim_method' in dsattrs and
                dsattrs['trim_method'].lower() == 'water level sl'):
            vel.attrs['note'] = ('Velocity bins trimmed if out of water or if '
                                 'side lobes intersect sea surface')

    def add_attributes(var, dsattrs):
        var.attrs.update(
            {'serial_number': dsattrs['AQDSerial_Number'],
             'initial_instrument_height': dsattrs['initial_instrument_height'],
             'nominal_instrument_depth': dsattrs['nominal_instrument_depth'],
             'height_depth_units': 'm',
             'sensor_type': dsattrs['INST_TYPE']})
        var.encoding['_FillValue'] = 1e35

    ds.attrs['COMPOSITE'] = np.int32(0)

    if utils.is_cf(ds):
        ds.attrs['featureType'] = 'timeSeriesProfile'

    # Update attributes for EPIC and STG compliance
    ds = utils.ds_coord_no_fillvalue(ds)

    ds['time'].attrs.update(
        {'standard_name': 'time',
         'axis': 'T'})

    if 'epic_time' in ds:
        ds['epic_time'].attrs.update(
            {'units': 'True Julian Day',
             'type': 'EVEN',
             'epic_code': 624})

    if 'epic_time2' in ds:
        ds['epic_time2'].attrs.update(
            {'units': 'msec since 0:00 GMT',
             'type': 'EVEN',
             'epic_code': 624})

    if 'epic_time_2d' in ds:
        ds['epic_time_2d'].attrs = ds['epic_time'].attrs
        ds['epic_time_2d'].attrs['type'] = 'UNEVEN'
    if 'epic_time2_2d' in ds:
        ds['epic_time2_2d'].attrs = ds['epic_time2'].attrs
        ds['epic_time2_2d'].attrs['type'] = 'UNEVEN'

    ds['depth'].attrs.update(
        {'units': 'm',
         'long_name': 'mean water depth',
         'initial_instrument_height': ds.attrs['initial_instrument_height'],
         'nominal_instrument_depth': ds.attrs['nominal_instrument_depth'],
         'epic_code': 3})

    if 'NAVD88_ref' in ds.attrs:
        ds['depth'].attrs['VERT_DATUM'] = 'NAVD88'
        ds['depth'].attrs['NOTE'] = ('Computed as platform depth [m NAVD88] '
                                     '- initial_instrument_height - bin '
                                     'distance from transducer')

    if not waves:
        ds['u_1205'].attrs.update(
            {'name': 'u',
             'long_name': 'Eastward Velocity',
             'generic_name': 'u',
             'epic_code': 1205})

        ds['v_1206'].attrs.update(
            {'name': 'v',
             'long_name': 'Northward Velocity',
             'generic_name': 'v',
             'epic_code': 1206})

        ds['w_1204'].attrs.update(
            {'name': 'w',
             'long_name': 'Vertical Velocity',
             'generic_name': 'w',
             'epic_code': 1204})

        ds['AGC_1202'].attrs.update(
            {'units': 'counts',
             'name': 'AGC',
             'long_name': 'Average Echo Intensity',
             'generic_name': 'AGC',
             'epic_code': 1202})

    elif waves:
        ds['vel1_1277'].attrs.update(
            {'units': 'mm/s',
             'long_name': 'Beam 1 Velocity',
             'generic_name': 'vel1',
             'epic_code': 1277})

        ds['vel2_1278'].attrs.update(
            {'units': 'mm/s',
             'long_name': 'Beam 2 Velocity',
             'generic_name': 'vel2',
             'epic_code': 1278})

        ds['vel3_1279'].attrs.update(
            {'units': 'mm/s',
             'long_name': 'Beam 3 Velocity',
             'generic_name': 'vel3',
             'epic_code': 1279})

        ds['AGC1_1221'].attrs.update(
            {'units': 'counts',
             'long_name': 'Echo Intensity (AGC) Beam 1',
             'generic_name': 'AGC1',
             'epic_code': 1221})

        ds['AGC2_1222'].attrs.update(
            {'units': 'counts',
             'long_name': 'Echo Intensity (AGC) Beam 2',
             'generic_name': 'AGC2',
             'epic_code': 1222})

        ds['AGC3_1223'].attrs.update(
            {'units': 'counts',
             'long_name': 'Echo Intensity (AGC) Beam 3',
             'generic_name': 'AGC3',
             'epic_code': 1223})

        ds.attrs['COORD_SYSTEM'] = 'GEOGRAPHIC + sample'

    ds['P_1'].attrs.update(
        {'units': 'dbar',
         'name': 'P',
         'long_name': 'Pressure',
         'generic_name': 'depth',
         'epic_code': 1})  # TODO: is this generic name correct?

    if 'P_1ac' in ds:
        ds['P_1ac'].attrs.update(
            {'units': 'dbar',
             'name': 'Pac',
             'long_name': 'Corrected pressure'})
        if 'P_1ac_note' in ds.attrs:
            ds['P_1ac'].attrs.update({'note': ds.attrs['P_1ac_note']})

        add_attributes(ds['P_1ac'], ds.attrs)

        ds.attrs['history'] = ('Atmospheric pressure compensated. ' +
                               ds.attrs['history'])

    ds['bin_depth'].attrs.update(
        {'units': 'm',
         'name': 'bin depth',
         'long_name': 'bin depth'})

    if 'P_1ac' in ds:
        if waves:
            ds['bin_depth'].attrs['note'] = ('Actual depth time series of '
                                             'wave burst bin depths. '
                                             'Calculated as corrected '
                                             'pressure (P_1ac) - bindist.')
        else:
            ds['bin_depth'].attrs['note'] = ('Actual depth time series of '
                                             'velocity bins. Calculated as '
                                             'corrected pressure (P_1ac) - '
                                             'bindist.')
    else:
        ds['bin_depth'].attrs.update(
            {'note': ('Actual depth time series of velocity bins. Calculated '
                      'as pressure (P_1) - bindist.')})

    ds['Tx_1211'].attrs.update(
        {'units': 'C',
         'name': 'Tx',
         'long_name': 'Instrument Transducer Temperature',
         'generic_name': 'temp',
         'epic_code': 1211})

    ds['Hdg_1215'].attrs.update(
        {'units': 'degrees',
         'name': 'Hdg',
         'long_name': 'Instrument Heading',
         'generic_name': 'hdg',
         'epic_code': 1215})

    if 'magnetic_variation_at_site' in ds.attrs:
        ds['Hdg_1215'].attrs['note'] = ('Heading is degrees true. Converted '
                                        'from magnetic with magnetic variation'
                                        ' of %f.' %
                                        ds.attrs['magnetic_variation_at_site'])
    elif 'magnetic_variation' in ds.attrs:
        ds['Hdg_1215'].attrs['note'] = ('Heading is degrees true. Converted '
                                        'from magnetic with magnetic variation'
                                        ' of %f.' %
                                        ds.attrs['magnetic_variation'])

    ds['Ptch_1216'].attrs.update(
        {'units': 'degrees',
         'name': 'Ptch',
         'long_name': 'Instrument Pitch',
         'generic_name': 'ptch',
         'epic_code': 1216})

    ds['Roll_1217'].attrs.update(
        {'units': 'degrees',
         'name': 'Roll',
         'long_name': 'Instrument Roll',
         'generic_name': 'roll',
         'epic_code': 1217})

    ds['Bat_106'].attrs.update({'units': 'V',
                                'long_name': 'Battery voltage',
                                'epic_code': 106})

    ds['bindist'].attrs.update(
        {'units': 'm',
         'long_name': 'distance from transducer head',
         'blanking_distance': ds.attrs['AQDBlankingDistance'],
         'note': ('distance is along profile from instrument '
                  'head to center of bin')})

    if not waves:
        for v in ['AGC_1202', 'u_1205', 'v_1206', 'w_1204']:
            add_attributes(ds[v], ds.attrs)
        for v in ['u_1205', 'v_1206', 'w_1204']:
            add_vel_attributes(ds[v], ds.attrs)
    elif waves:
        for v in ['vel1_1277',
                  'vel2_1278',
                  'vel3_1279',
                  'AGC1_1221',
                  'AGC2_1222',
                  'AGC3_1223']:
            add_attributes(ds[v], ds.attrs)

    for v in ['P_1',
              'Tx_1211',
              'Hdg_1215',
              'Ptch_1216',
              'Roll_1217',
              'Bat_106',
              'bin_depth',
              'bindist']:
        add_attributes(ds[v], ds.attrs)

    return ds
