from __future__ import division, print_function

import pandas as pd
import xarray as xr
import warnings
from ..core import utils

def prf_to_cdf(metadata):
    """Main load file"""

    # TODO: clock drift code
    # TODO: Move time to center of ensemble??
    # TODO: logmeta code
    # FIXME: Are fillvalues being inserted properly?

    basefile = metadata['basefile']

    # get instrument metadata from the HDR file
    instmeta = read_aqd_hdr(basefile)

    metadata['instmeta'] = instmeta

    print("Loading ASCII files")

    # Load sensor data
    ds = load_sen(basefile)

    # write out metadata first, then deal exclusively with xarray attrs
    ds = utils.write_metadata(ds, metadata)
    ds = utils.write_metadata(ds, metadata['instmeta'])

    del metadata
    del instmeta

    # Deal with metadata peculiarities
    ds = check_attrs(ds)

    ds = check_orientation(ds)

    # Load amplitude and velocity data
    ds = load_amp_vel(ds, basefile)

    # Compute time stamps
    ds = compute_time(ds)

    # configure file
    cdf_filename = ds.attrs['filename'] + '-raw.cdf'

    ds = update_attrs(ds)

    # need to drop datetime
    ds = ds.drop('datetime')

    ds.to_netcdf(cdf_filename, unlimited_dims='time')

    print('Finished writing data to %s' % cdf_filename)

    return ds


def load_sen(basefile):
    """Load data from .sen file"""

    senfile = basefile + '.sen'

    # read csv and parse dates
    # https://stackoverflow.com/questions/27112591/parsing-year-month-day-hour-minute-second-in-python
    def parse(year, month, day, hour, minute, second):
        return year + '-' + month + '-' + day + ' ' + hour + ':' + minute + ':' + second

    SEN = pd.read_csv(senfile, header=None, delim_whitespace=True,
        parse_dates={'datetime': [2, 0, 1, 3, 4, 5]}, date_parser=parse, usecols=[0, 1, 2, 3, 4, 5, 8, 10, 11, 12, 13, 14, 15, 16])

    # rename columns from numeric to human-readable
    SEN.rename(columns={10: 'Heading', 11: 'Pitch', 12: 'Roll', 13:'Pressure', 14:'Temperature', 8: 'Battery'}, inplace=True)

    # Look for analog data TODO
    SEN.rename(columns={15: 'AnalogInput1'}, inplace=True)
    SEN['AnalogInput1'] = SEN['AnalogInput1'] * 5 / 65535
    SEN.rename(columns={16: 'AnalogInput2'}, inplace=True)
    SEN['AnalogInput2'] = SEN['AnalogInput2'] * 5 / 65535

    # create xarray Dataset
    RAW = xr.Dataset.from_dataframe(SEN)
    RAW = RAW.rename({'index': 'time'})
    RAW['time'] = RAW['datetime']

    return RAW


def check_orientation(ds, waves=False):
    """Check instrument orientation and create variables that depend on this"""

    print('Insrument orientation:', ds.attrs['orientation'])
    print('Center_first_bin = %f' % ds.attrs['center_first_bin'])
    print('bin_size = %f' % ds.attrs['bin_size'])
    print('bin_count = %f' % ds.attrs['bin_count'])
    # TODO: these values are already in the HDR file...
    if not waves:
        bindist = np.linspace(ds.attrs['center_first_bin'],
                                     (ds.attrs['center_first_bin'] + ((ds.attrs['bin_count'] - 1) * ds.attrs['bin_size'])),
                                     num=ds.attrs['bin_count'])
    else:
        bindist = ds['cellpos'][0]

    if ds.attrs['orientation'] == 'UP':
        print('User instructed that instrument was pointing UP')
        # depth, or distance below surface, is a positive number below the
        # surface, negative above the surface, for CMG purposes and consistency with ADCP
        depth = (ds.attrs['WATER_DEPTH'] - ds.attrs['transducer_offset_from_bottom']) - bindist
        Depth_NOTE = 'user reports uplooking bin depths = water_depth - transducer offset from bottom - bindist' # TODO: this is never used
    elif ds.attrs['orientation'] == 'DOWN':
        print('User instructed that instrument was pointing DOWN')
        depth = (ds.attrs['WATER_DEPTH'] - ds.attrs['transducer_offset_from_bottom']) + bindist
        Depth_NOTE = 'user reports downlooking bin depths = water_depth - transducer_offset_from_bottom + bindist' # TODO: this is never used

    if not waves:
        ds['bindist'] = xr.DataArray(bindist, dims=('bindist'), name='bindist')
        ds['depth'] = xr.DataArray(depth, dims=('bindist'), name='depth')
    else:
        ds['bindist'] = xr.DataArray([bindist], dims=('bindist'), name='bindist')
        ds['depth'] = xr.DataArray([depth], dims=('bindist'), name='depth')

    return ds


def check_attrs(ds, waves=False):

    # Add some metadata originally in the run scripts
    ds.attrs['nominal_sensor_depth_note'] = 'WATER_DEPTH - initial_instrument_height'
    ds.attrs['nominal_sensor_depth'] = ds.attrs['WATER_DEPTH'] - ds.attrs['initial_instrument_height']
    ds.attrs['transducer_offset_from_bottom'] = ds.attrs['initial_instrument_height']

    # % now verify the global metadata for standard EPIC and cmg stuff
    # % everything in metadata and instmeta get written as global attributes
    # % these also get copied to the .nc file
    if 'initial_instrument_height' not in ds.attrs or np.isnan(ds.attrs['initial_instrument_height']):
        ds.attrs['initial_instrument_height'] = 0

    ds.attrs['serial_number'] = ds.attrs['AQDSerial_Number']

    # update metadata from Aquadopp header to CMG standard so that various
    # profilers have the same attribute wording.  Redundant, but necessary
    if not waves:
        ds.attrs['bin_count'] = ds.attrs['AQDNumberOfCells']
        ds.attrs['bin_size'] = ds.attrs['AQDCellSize'] / 100 # from cm to m
        ds.attrs['blanking_distance'] = ds.attrs['AQDBlankingDistance'] # already in m
        # Nortek lists the distance to the center of the first bin as the blanking
        # distance plus one cell size
        ds.attrs['center_first_bin'] = ds.attrs['blanking_distance'] + ds.attrs['bin_size'] # in m
    else:
        ds.attrs['bin_count'] = 1 # only 1 wave bin
        ds.attrs['bin_size'] = ds.attrs['WaveCellSize'] # already in m
        ds.attrs['blanking_distance'] = ds.attrs['AQDBlankingDistance'] # already in m
        # need to set center_first_bin after return in main calling function

    ds.attrs['salinity_set_by_user'] = ds.attrs['AQDSalinity']
    ds.attrs['salinity_set_by_user_units'] = 'ppt'

    ds.attrs['frequency'] = ds.attrs['AQDFrequency']
    ds.attrs['beam_width'] = ds.attrs['AQDBeamWidth']
    ds.attrs['beam_pattern'] = ds.attrs['AQDBeamPattern']
    ds.attrs['beam_angle'] = ds.attrs['AQDBeamAngle']
    # instmeta['AQDHeadRotation'] = metadata.pop('head_rotation') # also deletes this key. Not sure why the Matlab file does this, maybe need TODO look into this

    # TODO: figure these out
    # metadata['insterr'] = instmeta['error']
    # metadata['inststat'] = instmeta['status']
    # metadata['instorient'] = instmeta['orient']

    ds.attrs['INST_TYPE'] = 'Nortek Aquadopp Profiler';

    return ds


def update_attrs(ds, waves=False):
    """Define dimensions and variables in NetCDF file"""

    ds['lat'] = xr.DataArray([ds.attrs['latitude']], dims=('lat'), name='lat')
    ds['lon'] = xr.DataArray([ds.attrs['longitude']], dims=('lon'), name='lon')

    ds['TransMatrix'] = xr.DataArray(ds.attrs['AQDTransMatrix'], dims=('Tmatrix', 'Tmatrix'), name='TransMatrix')

    ds['time'].attrs.update({'standard_name': 'time',
        'axis': 'T'})

    ds['lat'].attrs.update({'units': 'degree_north',
        'long_name': 'Latitude',
        'epic_code': 500})

    ds['lon'].attrs.update({'units': 'degree_east',
        'long_name': 'Longitude',
        'epic_code': 502})

    ds['bindist'].attrs.update({'units': 'm',
        'long_name': 'distance from transducer head',
        'bin_size': ds.attrs['bin_size'],
        'center_first_bin': ds.attrs['center_first_bin'],
        'bin_count': ds.attrs['bin_count'],
        'transducer_offset_from_bottom': ds.attrs['transducer_offset_from_bottom']})

    ds['Temperature'].attrs.update({'units': 'C',
        'long_name': 'Temperature',
        'generic_name': 'temp'})

    ds['Pressure'].attrs.update({'units': 'dbar',
        'long_name': 'Pressure',
        'generic_name': 'press',
        'note': 'Raw pressure from instrument, not corrected for changes in atmospheric pressure'})

    for n in [1, 2, 3]:
        ds['VEL' + str(n)].attrs.update({'units': 'cm/s',
            'Type': 'scalar',
            'transducer_offset_from_bottom': ds.attrs['transducer_offset_from_bottom']})
        ds['AMP' + str(n)].attrs.update({'long_name': 'Beam ' + str(n) + ' Echo Amplitude',
            'units': 'counts',
            'Type': 'scalar',
            'transducer_offset_from_bottom': ds.attrs['transducer_offset_from_bottom'] })

    if not waves:
        veltxt = 'current velocity'
    else:
        veltxt = 'wave-burst velocity'

    if ds.attrs['AQDCoordinateSystem'] == 'ENU':
        ds['VEL1'].attrs.update({'long_name': 'Eastward ' + veltxt})
        ds['VEL2'].attrs.update({'long_name': 'Northward ' + veltxt})
        ds['VEL3'].attrs.update({'long_name': 'Vertical ' + veltxt})
    elif ds.attrs['AQDCoordinateSystem'] == 'XYZ':
        ds['VEL1'].attrs.update({'long_name': veltxt.capitalize() + ' in X Direction'})
        ds['VEL2'].attrs.update({'long_name': veltxt.capitalize() + ' in Y Direction'})
        ds['VEL3'].attrs.update({'long_name': veltxt.capitalize() + ' in Z Direction'})
    elif ds.attrs['AQDCoordinateSystem'] == 'BEAM':
        ds['VEL1'].attrs.update({'long_name': 'Beam 1 ' + veltxt})
        ds['VEL2'].attrs.update({'long_name': 'Beam 2 ' + veltxt})
        ds['VEL3'].attrs.update({'long_name': 'Beam 3 ' + veltxt})

    ds['Battery'].attrs.update({'units': 'Volts',
        'long_name': 'Battery Voltage'})

    ds['Pitch'].attrs.update({'units': 'degrees',
        'long_name': 'Instrument Pitch'})

    ds['Roll'].attrs.update({'units': 'degrees',
        'long_name': 'Instrument Roll'})

    ds['Heading'].attrs.update({'units': 'degrees',
        'long_name': 'Instrument Heading',
        'datum': 'magnetic north'})

    ds['depth'].attrs.update({'units': 'm',
        'long_name': 'mean water depth',
        'bin_size': ds.attrs['bin_size'],
        'center_first_bin': ds.attrs['center_first_bin'],
        'bin_count': ds.attrs['bin_count'],
        'transducer_offset_from_bottom': ds.attrs['transducer_offset_from_bottom']})

    ds['TransMatrix'].attrs.update({'long_name': 'Transformation Matrix for this Aquadopp'})

    return ds

    # RAW['AnalogInput1']

    # with Dataset(cdf_filename, 'w', format='NETCDF4', clobber=True) as rg:
    #
    #     # write out EPIC metadata
    #     write_metadata(rg, RAW['instmeta']) #TODO
    #
    #     if waves:
    #         burstid = rg.createVariable('burst', 'i', ('time',), fill_value=False)
    #         burstid.units = 'count'
    #         burstid.long_name = 'Record Number'
    #
    #
    #
    #     if waves:
    #         AMP1id = rg.createVariable('AMP1', 'f', ('sample', 'time',), fill_value=False)
    #     else:
    #         AMP1id = rg.createVariable('AMP1', 'f', ('depth', 'time',), fill_value=False)
    #
    # FIXME: add analoginput stuff
    #     for n in ['1', '2']:
    #         if 'AnalogInput' + n in metadata:
    #             Anaid = rg.createVariable('AnalogInput' + n, 'f', ('time',), fill_value=False)
    #             Anaid.units = 'Volts'
    #             Anaid.sensor_type = metadata['AnalogInput1']['sensor_type']
    #             Anaid.sensor_manufacturer = metadata['AnalogInput1']['sensor_manufacturer']
    #             Anaid.sensor_model = metadata['AnalogInput1']['sensor_model']
    #             Anaid.serial_number = metadata['AnalogInput1']['serial_number']
    #
    #             if 'initial_sensor_height' in metadata['AnalogInput' + n]:
    #                 # TODO
    #                 # metadata.AnalogInput1.nominal_sensor_depth = metadata.WATER_DEPTH - metadata.AnalogInput1.initial_sensor_height;
    #                 # netcdf.putAtt(ncid,Ana1id,'initial_sensor_height',metadata.AnalogInput1.initial_sensor_height);
    #                 # netcdf.putAtt(ncid,Ana1id,'nominal_sensor_depth',metadata.AnalogInput1.nominal_sensor_depth);
    #                 continue
    #             elif 'nominal_sensor_depth' in metadata['AnalogInput' + n]: # TODO: should be another if not elif??
    #                 # netcdf.putAtt(ncid,Ana1id,'nominal_sensor_depth',metadata.AnalogInput1.nominal_sensor_depth);
    #                 # metadata.AnalogInput1.initial_sensor_height = metadata.WATER_DEPTH - metadata.AnalogInput1.nominal_sensor_depth;
    #                 # netcdf.putAtt(ncid,Ana1id,'initial_sensor_height',metadata.AnalogInput1.initial_sensor_height);
    #                 continue
    #
    #         # if isfield(metadata.AnalogInput1,'range'),
    #         #         netcdf.putAtt(ncid,Ana1id,'range',metadata.AnalogInput1.range);
    #         # end
    #         # if isfield(metadata.AnalogInput1.cals,'NTUcoef'),
    #         #         netcdf.putAtt(ncid,Ana1id,'NTUcoef',metadata.AnalogInput1.cals.NTUcoef);
    #         # end
    #         # if isfield(metadata.AnalogInput1.cals,'SEDcoef'),
    #         #         netcdf.putAtt(ncid,Ana1id,'SEDcoef',metadata.AnalogInput1.cals.SEDcoef);
    #         # end


def compute_time(ds, waves=False):
    """Compute Julian date and then time and time2 for use in netCDF file"""

    # shift times to center of ensemble
    if not waves:
        timeshift = ds.attrs['AQDAverageInterval']/2
    else:
        fs = float(ds.attrs['WaveSampleRate'].split()[0])
        timeshift = ds.attrs['WaveNumberOfSamples']/fs/2

    if timeshift.is_integer():
        ds['time'] = ds['time'] + np.timedelta64(int(timeshift), 's')
        print('Time shifted by:', int(timeshift), 's')
    else:
        warnings.warn('time NOT shifted because not a whole number of seconds: %f s ***' % timeshift)

    # create Julian date
    ds['jd'] = ds['time'].to_dataframe().index.to_julian_date() + 0.5

    ds['epic_time'] = np.floor(ds['jd'])
    if np.all(np.mod(ds['epic_time'], 1) == 0): # make sure they are all integers, and then cast as such
        ds['epic_time'] = ds['epic_time'].astype(np.int32)
    else:
        warnings.warn('not all EPIC time values are integers; this will cause problems with time and time2')

    # TODO: Hopefully this is correct... roundoff errors on big numbers...
    ds['epic_time2'] = np.round((ds['jd'] - np.floor(ds['jd']))*86400000).astype(np.int32)

    return ds


def load_amp_vel(RAW, basefile):
    """Load amplitude and velocity data from the .aN and .vN files"""

    for n in [1, 2, 3]:
        afile = basefile + '.a' + str(n)
        RAW['AMP' + str(n)] = xr.DataArray(pd.read_csv(afile, header=None, delim_whitespace=True),
            dims=('time', 'bindist'), coords=[RAW['time'], RAW['bindist']])
        # RAW['AMP' + str(n)] = RAW['AMP' + str(n)].rename({'dim_0': 'time'})
        vfile = basefile + '.v' + str(n)
        v = pd.read_csv(vfile, header=None, delim_whitespace=True)
        # convert to cm/s
        RAW['VEL' + str(n)] = xr.DataArray(v * 100, dims=('time', 'bindist'), coords=[RAW['time'], RAW['bindist']])

    return RAW


def read_aqd_hdr(basefile):
    """
    Get instrument metadata from .hdr file
    Was formerly readAQDprfHeader.m
    """
    #
    # % TODO read and save all of the instrument settings in the .hdr file
    # % replacing strncmp with strfind, strncmp need the exact length of string,
    # % many of those were wrong and lots of metadata was going missing.
    hdrFile = basefile + '.hdr'
    f = open(hdrFile, 'r')
    row = ''

    Instmeta = {}

    while 'Hardware configuration' not in row:
        row = f.readline().rstrip()
        if 'Profile interval' in row:
            idx = row.find(' sec')
            Instmeta['AQDProfileInterval'] = float(row[38:idx])
        elif 'Number of cells' in row:
            Instmeta['AQDNumberOfCells'] = float(row[38:])
        elif row.find('Cell size', 0, 9) != -1: # required here to differentiate from the wave cell size
            idx = row.find(' cm')
            Instmeta['AQDCellSize'] = float(row[38:idx])
        elif 'Average interval' in row:
            idx = row.find(' sec')
            Instmeta['AQDAverageInterval'] = float(row[38:idx])
        elif 'Measurement load' in row:
            idx = row.find(' %')
            Instmeta['AQDMeasurementLoad'] = float(row[38:idx])
        elif 'Transmit pulse length' in row:
            idx = row.find(' m')
            Instmeta['AQDTransmitPulseLength'] = float(row[38:idx])
        elif 'Blanking distance' in row:
            idx = row.find(' m')
            Instmeta['AQDBlankingDistance'] = float(row[38:idx])
        elif 'Compass update rate' in row:
            idx = row.find(' sec')
            Instmeta['AQDCompassUpdateRate'] = float(row[38:idx])
        elif 'Wave measurements' in row:
            Instmeta['WaveMeasurements'] = row[38:]
        elif 'Wave - Powerlevel' in row:
            Instmeta['WavePower'] = row[38:]
        elif 'Wave - Interval' in row:
            idx = row.find(' sec')
            Instmeta['WaveInterval'] = float(row[38:idx])
        elif 'Wave - Number of samples' in row:
            Instmeta['WaveNumberOfSamples'] = float(row[38:])
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
        elif 'Powerlevel' in row: # TODO: WRONG, this is not analog powerlevel
            Instmeta['AQDAnalogPowerLevel'] = row[38:]
        elif 'Coordinate system' in row:
            Instmeta['AQDCoordinateSystem'] = row[38:]
        elif 'Sound speed' in row:
            Instmeta['AQDSoundSpeed'] = row[38:]
        elif 'Salinity' in row:
            Instmeta['AQDSalinity'] = row[38:]
        elif 'Number of beams' in row:
            Instmeta['AQDNumberOfBeams'] = float(row[38:])
        elif 'Number of pings per burst' in row:
            Instmeta['AQDNumberOfPingsPerBurst'] = float(row[38:])
        elif 'Software version' in row:
            Instmeta['AQDSoftwareVersion'] = row[38:] # can't be float, was wrong in m-file
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
    # while 'Data file format' not in row:
        row = f.readline().rstrip()
        if 'Pressure sensor' in row:
            Instmeta['AQDPressureSensor'] = row[38:]
        elif 'Compass' in row:
            Instmeta['AQDCompass'] = row[38:]
        elif 'Tilt sensor' in row:
            Instmeta['AQDTilt'] = row[38:]
        elif 'Head frequency' in row:
            idx = row.find(' kHz')
            Instmeta['AQDFrequency'] = float(row[38:idx])
        elif 'Number of beams' in row:
            Instmeta['AQDNumBeams'] = float(row[38:])
        elif 'Serial number' in row:
            Instmeta['AQDHeadSerialNumber'] = row[38:]
        elif 'Transformation matrix' in row:
            Instmeta['AQDTransMatrix'] = np.zeros((3,3))
            for n in np.arange(3):
                Instmeta['AQDTransMatrix'][n,:] = [float(x) for x in row[38:].split()]
                row = f.readline().rstrip()
        elif 'Pressure sensor calibration' in row:
            Instmeta['AQDPressureCal'] = row[38:]

    # % infer some things based on the Aquadopp brochure
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
    Instmeta['AQDBeamAngle'] = 25;
# Instmeta.AQDVelRange = 1000; % cm/s
# Instmeta.AQDTempRange = [-4 40];
# Instmeta.AQDPressRange = [0 100];
# % no tilt range given in AQD docs
#
# fclose(hdr);
#
# % %if waves data were not collected remove wave parameters from metadata
# % if strfind(Instmeta.AQDWaveStatus,'DISABLED',7)
# %     fields = {'AQDWavePower','AQDWaveInterval','AQDWaveSampleRate','AQDWaveNumberOfSamples'};
# %     Instmeta = rmfield(Instmeta,fields);
# % else
# % end

    return Instmeta


def insert_fill_values(RAW):
    """Insert fill values for nans"""

    print("Inserting fill values")
    for k in RAW:
        if k not in ['instmeta', 'time', 'time2', 'datetime'] and np.max(np.shape(RAW[k])) == np.max(np.shape(RAW['jd'])):
            nanind = np.where(np.isnan(RAW[k]))
            RAW[k][nanind] = 1e35

    return RAW
