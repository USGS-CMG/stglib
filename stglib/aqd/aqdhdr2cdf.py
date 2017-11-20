#!/usr/bin/env python

from __future__ import division, print_function
import numpy as np
import sys
sys.path.insert(0, '/Users/dnowacki/Documents/aqdlib')
import aqdlib
import aqdlib.qaqc
import pandas as pd
import xarray as xr
import warnings
import platform
import netCDF4
import inspect
import os


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
    RAW = load_sen(metadata)

    # Deal with metadata peculiarities
    metadata = check_metadata(metadata)

    RAW = check_orientation(RAW, metadata)

    # Load amplitude and velocity data
    RAW = load_amp_vel(RAW, basefile)

    # Compute time stamps
    RAW = compute_time(RAW, metadata)

    # configure file
    cdf_filename = metadata['filename'] + '-raw.cdf'

    # write out metadata
    RAW = write_metadata(RAW, metadata)
    RAW = write_metadata(RAW, metadata['instmeta'])

    update_attrs(cdf_filename, RAW, metadata)

    # need to drop datetime
    RAW = RAW.drop('datetime')

    RAW.to_netcdf(cdf_filename, unlimited_dims='time')

    print('Finished writing data to %s' % cdf_filename)

    return RAW


def load_sen(metadata):
    """Load data from .sen file"""

    senfile = metadata['basefile'] + '.sen'

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


def check_orientation(RAW, metadata, waves=False):
    """Check instrument orientation and create variables that depend on this"""

    print('Insrument orientation:', metadata['orientation'])
    print('Center_first_bin = %f' % metadata['center_first_bin'])
    print('bin_size = %f' % metadata['bin_size'])
    print('bin_count = %f' % metadata['bin_count'])
    # TODO: these values are already in the HDR file...
    if not waves:
        bindist = np.linspace(metadata['center_first_bin'],
                                     (metadata['center_first_bin'] + ((metadata['bin_count'] - 1) * metadata['bin_size'])),
                                     num=metadata['bin_count'])
    else:
        bindist = RAW['cellpos'][0]

    if metadata['orientation'] == 'UP':
        print('User instructed that instrument was pointing UP')
        # depth, or distance below surface, is a positive number below the
        # surface, negative above the surface, for CMG purposes and consistency with ADCP
        depth = (metadata['WATER_DEPTH'] - metadata['transducer_offset_from_bottom']) - bindist
        Depth_NOTE = 'user reports uplooking bin depths = water_depth - transducer offset from bottom - bindist' # TODO: this is never used
    elif metadata['orientation'] == 'DOWN':
        print('User instructed that instrument was pointing DOWN')
        depth = (metadata['WATER_DEPTH'] - metadata['transducer_offset_from_bottom']) + bindist
        Depth_NOTE = 'user reports downlooking bin depths = water_depth - transducer_offset_from_bottom + bindist' # TODO: this is never used

    if not waves:
        RAW['bindist'] = xr.DataArray(bindist, dims=('bindist'), name='bindist')
        RAW['depth'] = xr.DataArray(depth, dims=('bindist'), name='depth')
    else:
        RAW['bindist'] = xr.DataArray([bindist], dims=('bindist'), name='bindist')
        RAW['depth'] = xr.DataArray([depth], dims=('bindist'), name='depth')

    return RAW


def check_metadata(metadata, waves=False):

    # Add some metadata originally in the run scripts
    metadata['nominal_sensor_depth_note'] = 'WATER_DEPTH - initial_instrument_height'
    metadata['nominal_sensor_depth'] = metadata['WATER_DEPTH'] - metadata['initial_instrument_height']
    metadata['transducer_offset_from_bottom'] = metadata['initial_instrument_height']

    # % now verify the global metadata for standard EPIC and cmg stuff
    # % everything in metadata and instmeta get written as global attributes
    # % these also get copied to the .nc file
    if 'initial_instrument_height' not in metadata or np.isnan(metadata['initial_instrument_height']):
        metadata['initial_instrument_height'] = 0

    metadata['serial_number'] = metadata['instmeta']['AQDSerial_Number']

    # update metadata from Aquadopp header to CMG standard so that various
    # profilers have the same attribute wording.  Redundant, but necessary
    if not waves:
        metadata['bin_count'] = metadata['instmeta']['AQDNumberOfCells']
        metadata['bin_size'] = metadata['instmeta']['AQDCellSize'] / 100 # from cm to m
        metadata['blanking_distance'] = metadata['instmeta']['AQDBlankingDistance'] # already in m
        # Nortek lists the distance to the center of the first bin as the blanking
        # distance plus one cell size
        metadata['center_first_bin'] = metadata['blanking_distance'] + metadata['bin_size'] # in m
    else:
        metadata['bin_count'] = 1 # only 1 wave bin
        metadata['bin_size'] = metadata['instmeta']['WaveCellSize'] # already in m
        metadata['blanking_distance'] = metadata['instmeta']['AQDBlankingDistance'] # already in m
        # need to set center_first_bin after return in main calling function

    metadata['salinity_set_by_user'] = metadata['instmeta']['AQDSalinity']
    metadata['salinity_set_by_user_units'] = 'ppt'

    metadata['frequency'] = metadata['instmeta']['AQDFrequency']
    metadata['beam_width'] = metadata['instmeta']['AQDBeamWidth']
    metadata['beam_pattern'] = metadata['instmeta']['AQDBeamPattern']
    metadata['beam_angle'] = metadata['instmeta']['AQDBeamAngle']
    # instmeta['AQDHeadRotation'] = metadata.pop('head_rotation') # also deletes this key. Not sure why the Matlab file does this, maybe need TODO look into this

    # TODO: figure these out
    # metadata['insterr'] = instmeta['error']
    # metadata['inststat'] = instmeta['status']
    # metadata['instorient'] = instmeta['orient']

    metadata['INST_TYPE'] = 'Nortek Aquadopp Profiler';

    return metadata


def update_attrs(cdf_filename, RAW, metadata, waves=False):
    """Define dimensions and variables in NetCDF file"""

    RAW['lat'] = xr.DataArray([metadata['latitude']], dims=('lat'), name='lat')
    RAW['lon'] = xr.DataArray([metadata['longitude']], dims=('lon'), name='lon')

    RAW['TransMatrix'] = xr.DataArray(metadata['instmeta']['AQDTransMatrix'], dims=('Tmatrix', 'Tmatrix'), name='TransMatrix')

    RAW['time'].attrs.update({'standard_name': 'time',
        'axis': 'T'})

    RAW['lat'].attrs.update({'units': 'degree_north',
        'long_name': 'Latitude',
        'epic_code': 500})

    RAW['lon'].attrs.update({'units': 'degree_east',
        'long_name': 'Longitude',
        'epic_code': 502})

    RAW['bindist'].attrs.update({'units': 'm',
        'long_name': 'distance from transducer head',
        'bin_size': metadata['bin_size'],
        'center_first_bin': metadata['center_first_bin'],
        'bin_count': metadata['bin_count'],
        'transducer_offset_from_bottom': metadata['transducer_offset_from_bottom']})

    RAW['Temperature'].attrs.update({'units': 'C',
        'long_name': 'Temperature',
        'generic_name': 'temp'})

    RAW['Pressure'].attrs.update({'units': 'dbar',
        'long_name': 'Pressure',
        'generic_name': 'press',
        'note': 'raw pressure from instrument, not corrected for changes in atmospheric pressure'})

    for n in [1, 2, 3]:
        RAW['VEL' + str(n)].attrs.update({'units': 'cm/s',
            'Type': 'scalar',
            'transducer_offset_from_bottom': metadata['transducer_offset_from_bottom']})
        RAW['AMP' + str(n)].attrs.update({'long_name': 'Beam ' + str(n) + ' Echo Amplitude',
            'units': 'counts',
            'Type': 'scalar',
            'transducer_offset_from_bottom': metadata['transducer_offset_from_bottom'] })

    if not waves:
        veltxt = 'current velocity'
    else:
        veltxt = 'wave-burst velocity'

    if metadata['instmeta']['AQDCoordinateSystem'] == 'ENU':
        RAW['VEL1'].attrs.update({'long_name': 'Eastward ' + veltxt})
        RAW['VEL2'].attrs.update({'long_name': 'Northward ' + veltxt})
        RAW['VEL3'].attrs.update({'long_name': 'Vertical ' + veltxt})
    elif metadata['instmeta']['AQDCoordinateSystem'] == 'XYZ':
        RAW['VEL1'].attrs.update({'long_name': veltxt.capitalize() + ' in X Direction'})
        RAW['VEL2'].attrs.update({'long_name': veltxt.capitalize() + ' in Y Direction'})
        RAW['VEL3'].attrs.update({'long_name': veltxt.capitalize() + ' in Z Direction'})
    elif metadata['instmeta']['AQDCoordinateSystem'] == 'BEAM':
        RAW['VEL1'].attrs.update({'long_name': 'Beam 1 ' + veltxt})
        RAW['VEL2'].attrs.update({'long_name': 'Beam 2 ' + veltxt})
        RAW['VEL3'].attrs.update({'long_name': 'Beam 3 ' + veltxt})

    RAW['Battery'].attrs.update({'units': 'Volts',
        'long_name': 'Battery Voltage'})

    RAW['Pitch'].attrs.update({'units': 'degrees',
        'long_name': 'Instrument Pitch'})

    RAW['Roll'].attrs.update({'units': 'degrees',
        'long_name': 'Instrument Roll'})

    RAW['Heading'].attrs.update({'units': 'degrees',
        'long_name': 'Instrument Heading',
        'datum': 'magnetic north'})

    RAW['depth'].attrs.update({'units': 'm',
        'long_name': 'mean water depth',
        'bin_size': metadata['bin_size'],
        'center_first_bin': metadata['center_first_bin'],
        'bin_count': metadata['bin_count'],
        'transducer_offset_from_bottom': metadata['transducer_offset_from_bottom']})

    RAW['TransMatrix'].attrs.update({'long_name': 'Transformation Matrix for this Aquadopp'})

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


def write_metadata(ds, metadata):
    """Write out all metadata to CDF file"""

    for k in metadata:
        if k != 'instmeta': # don't want to write out instmeta dict, call it separately
            ds.attrs.update({k: metadata[k]})

    f = os.path.basename(inspect.stack()[1][1])

    ds.attrs.update({'history': 'Processed using ' + f + ' with Python ' +
        platform.python_version() + ', xarray ' + xr.__version__ + ', NumPy ' +
        np.__version__ + ', netCDF4 ' + netCDF4.__version__})

    return ds


def compute_time(RAW, metadata, waves=False):
    """Compute Julian date and then time and time2 for use in netCDF file"""

    # shift times to center of ensemble
    if not waves:
        timeshift = metadata['instmeta']['AQDAverageInterval']/2
    else:
        fs = float(metadata['instmeta']['WaveSampleRate'].split()[0])
        timeshift = metadata['instmeta']['WaveNumberOfSamples']/fs/2

    if timeshift.is_integer():
        RAW['time'] = RAW['time'] + np.timedelta64(int(timeshift), 's')
        print('Time shifted by:', int(timeshift), 's')
    else:
        warnings.warn('time NOT shifted because not a whole number of seconds: %f s ***' % timeshift)

    # create Julian date
    RAW['jd'] = RAW['time'].to_dataframe().index.to_julian_date() + 0.5

    RAW['epic_time'] = np.floor(RAW['jd'])
    if np.all(np.mod(RAW['epic_time'], 1) == 0): # make sure they are all integers, and then cast as such
        RAW['epic_time'] = RAW['epic_time'].astype(np.int32)
    else:
        warnings.warn('not all EPIC time values are integers; this will cause problems with time and time2')

    # TODO: Hopefully this is correct... roundoff errors on big numbers...
    RAW['epic_time2'] = np.round((RAW['jd'] - np.floor(RAW['jd']))*86400000).astype(np.int32)

    return RAW


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
            RAW[k][nanind] = aqdlib.DOUBLE_FILL

    return RAW


def main():
    import sys
    sys.path.insert(0, '/Users/dnowacki/Documents/aqdlib')
    import aqdlib
    import argparse
    import yaml

    parser = argparse.ArgumentParser(description='Convert Aquadopp text files to raw .cdf format. Run this script from the directory containing Aquadopp files')
    # parser.add_argument('basename', help='base name (without extension) of the Aquadopp text files')
    parser.add_argument('gatts', help='path to global attributes file (gatts formatted)')
    parser.add_argument('config', help='path to ancillary config file (YAML formatted)')
    # parser.add_argument('-v', '--verbose', help='increase output verbosity', action='store_true')

    args = parser.parse_args()

    # initialize metadata from the globalatts file
    metadata = aqdlib.read_globalatts(args.gatts)

    # Add additional metadata from metadata config file
    config = yaml.safe_load(open(args.config))

    for k in config:
        metadata[k] = config[k]

    RAW = prf_to_cdf(metadata)

if __name__ == '__main__':
    main()
