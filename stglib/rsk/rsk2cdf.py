from __future__ import division, print_function
import sqlite3
import numpy as np
import xarray as xr
import pandas as pd
from ..core import utils


def rsk_to_cdf(metadata):
    """
    Main function to load data from RSK file and save to raw .CDF
    """

    ds = rsk_to_xr(metadata)

    print("Writing to raw netCDF")

    ds.to_netcdf(ds.attrs['filename'] + '-raw.cdf')

    print("Done")

    return ds


def init_connection(rskfile):
    """Initialize an sqlite3 connection and return a cursor"""

    conn = sqlite3.connect(rskfile)
    return conn.cursor()


def rsk_to_xr(metadata):
    """
    Load data from RSK file and generate an xarray Dataset
    """

    rskfile = metadata.pop('basefile') + '.rsk'

    ds = xr.Dataset()

    ds = utils.write_metadata(ds, metadata)

    print(('Loading from sqlite file %s; '
           'this may take a while for large datasets') % rskfile)

    conn = init_connection(rskfile)
    
    # note that this code only reads the first deployment in an .rsk file

    # the field names differ for different RBRs, to we need to know
    # if this is a Virtuoso Tu because the data is very different from
    # the DWave for which this was written
    if 'instrument_type' in metadata and metadata['instrument_type'] == "Virtuoso Tu":
    
        ########  table: deployments  ########
        #['deploymentID', 'instrumentID', 'comment', 'loggerStatus', 'loggerTimeDrift', 'timeOfDownload', 'name', 'sampleSize']
        #(1, 1, '4th attempt\r\nNow with Ruskin v2.3.1', 'stopped', -582653584941, 946684866123, '\\054276_20180618_1614.rsk', 156864)

        # right now we are assuming a database with only on deployment.
        # if there are more than one deployment, only the first is read
        
        # some important keys into the RBR organization of their database
        # we do not use these right now - and they may be important if there is 
        # more than one deployment in an rsk file
        deploymentID = int(conn.execute(
            "select deploymentID from deployments").fetchall()[0][0])
        ds.attrs['RBR_deployment_ID'] = deploymentID
        instrumentID = int(conn.execute(
            "select instrumentID from deployments").fetchall()[0][0])
        ds.attrs['RBR_instrument_ID'] = instrumentID
        nrecords = int(conn.execute(
            "select samplesize from deployments").fetchall()[0][0])
        ds.attrs['RBR_number_of_records'] = nrecords
        ds.attrs['RBR_loggerTimeDrift'] = int(conn.execute(
            "select loggerTimeDrift from deployments").fetchall()[0][0])
        ds.attrs['comment'] = str(conn.execute(
            "select comment from deployments").fetchall()[0][0])

        conn.execute("SELECT * FROM data")
        data = conn.fetchall()
        print("Done fetching data")
        d = np.asarray(data)

        # Get record interval, which is continuous mode defined in schedules
        sampling_mode = str(conn.execute(
            "select mode from schedules").fetchall()[0][0])
        if "continuous" in sampling_mode:
            ds.attrs['samples_per_burst'] = 1
            # for burst data
            # samplingcount = ds.attrs['samples_per_burst']
        ds.attrs['sample_interval'] = int(conn.execute(
            "select samplingPeriod from "+sampling_mode).fetchall()[0][0]) / 1000
        #ds.attrs['burst_interval'] = ds.attrs['sample_interval']
        
        ########  table: instruments  ########
        #['instrumentID', 'serialID', 'model', 'firmwareVersion', 'firmwareType']
        #(1, 54276, 'RBRvirtuoso', '1.41', 103)
        ds.attrs['serial_number'] = str(conn.execute(
            "select serialID from instruments").fetchall()[0][0])
        ds.attrs['INST_TYPE'] = str(conn.execute(
            "select model from instruments").fetchall()[0][0])
        ds.attrs['RBR_firmwareVersion'] = str(conn.execute(
            "select firmwareVersion from instruments").fetchall()[0][0])
        ds.attrs['RBR_firmwareType'] = str(conn.execute(
            "select firmwareType from instruments").fetchall()[0][0])
        ds.attrs['RBR_ruskinVersion'] = str(conn.execute(
            "select ruskinVersion from appSettings").fetchall()[0][0])
        ########  table: ranging  ########
        #['instrumentID', 'channelID', 'channelOrder', 'mode', 'gain', 'availableGains']
        #(1, 1, 1, 'Auto', 1.0, '1.0|5.0|20.0|100.0')
        ds.attrs['RBR_ranging_mode'] = str(conn.execute(
            "select mode from ranging").fetchall()[0][0])
        ds.attrs['RBR_ranging_gain'] = str(conn.execute(
            "select gain from ranging").fetchall()[0][0])
        ds.attrs['RBR_ranging_availableGains'] = str(conn.execute(
            "select availableGains from ranging").fetchall()[0][0])
        ########  table: coefficients  ########
        #['calibrationID', 'key', 'value']
        #(1, 'c0', '3483.989')
        #(1, 'c1', '-4756.7354')
        coeffs = conn.execute("select key from coefficients").fetchall()
        coeff_vals = conn.execute("select value from coefficients").fetchall()
        for idx in range(len(coeffs)):
            ds.attrs['RBR_coefficients_'+coeffs[idx][0]] = coeff_vals[idx][0]
            # Some other specifics to be aware of and where they are in the rsk schema
        ########  table: channels  ########
        #['channelID', 'shortName', 'longName', 'units', 'longNamePlainText', 'unitsPlainText', 'isMeasured', 'isDerived']
        #(1, 'turb00', 'Turbidity', 'NTU', 'Turbidity', 'NTU', 1, 0)

        # what variable is this going to be in the netcdf?
        # Let's not hard code it so we can change it easily later
        # for some other sensor on some other single channel RBR
        varname = 'Tu'
        epic_varname = 'Turb'
        epic_name = 'Turbidity'
        
        # d is a matrix shaped as [time, measurement type per column]
        # simplified case of burst data - still reshape
        a = {}
        a['unixtime'] = d[:, 0].copy()
        a[varname] = d[:, 1].copy()
        # sort by time (not sorted for some reason)
        sort = np.argsort(a['unixtime'])
        a['unixtime'] = a['unixtime'][sort]
        a[varname] = a[varname][sort]
        
        # for continuous data
        times = pd.to_datetime(a['unixtime'][:], unit='ms')

        depth = np.array(ds.attrs['WATER_DEPTH']-ds.attrs['initial_instrument_height'])
        lat = ds.attrs['latitude']
        lon = ds.attrs['longitude']
    
        print('we are up to writing Turb')

        # version for continuous data set
        ds[epic_varname] = xr.DataArray(
            a[varname][:],
            coords=[times],
            dims=('time'),
            name=epic_name,
            attrs={
                'long_name': str(conn.execute(
                    "select LongName from channels").fetchall()[0][0]),
                'units': str(conn.execute(
                    "select units from channels").fetchall()[0][0]),
                'initial_instrument_height': ds.attrs['initial_instrument_height'],
                'serial_number': ds.attrs['serial_number']},
            encoding={'_FillValue': 1e35})
            
        ds['depth'] = xr.DataArray(
            [depth],
            dims=('depth'),
            name='depth',
            attrs={'units': 'm',
                   'type': 'EVEN',
                   'epic_code': 3})
        
    else:  # we presume this is a DWave
        conn.execute("SELECT tstamp, channel01 FROM burstdata")
        data = conn.fetchall()
        print("Done fetching data")
        d = np.asarray(data)
        # Get samples per burst
        try:
            # this seems to be used on older-style databases;
            # throws error on newer files
            samplingcount = conn.execute(
                "select samplingcount from schedules").fetchall()[0][0]
        except sqlite3.OperationalError:
            samplingcount = conn.execute(
                "select samplingcount from wave").fetchall()[0][0]
        ds.attrs['samples_per_burst'] = samplingcount

        try:
            samplingperiod = conn.execute(
                "select samplingperiod from schedules").fetchall()[0][0]
        except sqlite3.OperationalError:
            samplingperiod = conn.execute(
                "select samplingperiod from wave").fetchall()[0][0]
        ds.attrs['sample_interval'] = samplingperiod / 1000

        try:
            repetitionperiod = conn.execute(
                "select repetitionperiod from schedules").fetchall()[0][0]
        except sqlite3.OperationalError:
            repetitionperiod = conn.execute(
                "select repetitionperiod from wave").fetchall()[0][0]
        ds.attrs['burst_interval'] = repetitionperiod / 1000

        ds.attrs['burst_length'] = ds.attrs['samples_per_burst'] * \
            ds.attrs['sample_interval']
        ds.attrs['serial_number'] = str(conn.execute(
            "select serialID from instruments").fetchall()[0][0])
        ds.attrs['INST_TYPE'] = 'RBR Virtuoso d|wave'
    
        #conn.close()

        a = {}
        a['unixtime'] = d[:, 0].copy()
        a['pres'] = d[:, 1].copy()
        # sort by time (not sorted for some reason)
        sort = np.argsort(a['unixtime'])
        a['unixtime'] = a['unixtime'][sort]
        a['pres'] = a['pres'][sort]

        # get indices that end at the end of the final burst
        datlength = a['unixtime'].shape[0] - a['unixtime'].shape[0] % samplingcount

        # reshape
        for k in a:
            a[k] = a[k][:datlength].reshape(
                (int(datlength/samplingcount), samplingcount)
            )

        times = pd.to_datetime(a['unixtime'][:, 0], unit='ms')
        samples = np.arange(samplingcount)

        ds['P_1'] = xr.DataArray(
            a['pres'],
            coords=[times, samples],
            dims=('time', 'sample'),
            name='Pressure',
            attrs={
                'long_name': 'Pressure',
                'units': 'dbar',
                'epic_code': 1,
                'height_depth_units': 'm',
                'initial_instrument_height': ds.attrs['initial_instrument_height'],
                'serial_number': ds.attrs['serial_number']},
            encoding={'_FillValue': 1e35})

        ds['time'] = xr.DataArray(times, dims=('time'), name='time')

        ds['sample'] = xr.DataArray(samples, dims=('sample'), name='sample')

        ds['burst'] = xr.DataArray(np.arange(len(times)), dims=('time'),
                                   attrs={'long_name': 'burst number'})

    ds['lat'] = xr.DataArray(
        [ds.attrs['latitude']],
        dims=('lat'),
        name='lat',
        attrs={'units': 'degree_north',
               'long_name': 'Latitude',
               'epic_code': 500})

    ds['lon'] = xr.DataArray(
        [ds.attrs['longitude']],
        dims=('lon'),
        name='lon',
        attrs={'units': 'degree_east',
               'long_name': 'Longitude',
               'epic_code': 502})

    # need to add  time attrs after DataArrays have been combined into Dataset
    ds['time'].attrs.update({'standard_name': 'time', 'axis': 'T'})
    
    conn.close()

    return ds

    # # TODO: add the following??
    # # {'positive','down';
    # #                'long_name', 'Depth';
    # #                'axis','z';
    # #                'units', 'm';
    # #                'epic_code', 3};
    #
    # Pressid = rg.createVariable('Pressure', 'f', ('time','sample',),
    # Pressid.units = 'dbar'
    # Pressid.long_name = 'Pressure (dbar)'
    # Pressid.generic_name = 'press'
    # Pressid.note = 'raw pressure from instrument, not corrected for...
    # Pressid.epic_code = 1
    # Pressid.height_depth_units = 'm'
