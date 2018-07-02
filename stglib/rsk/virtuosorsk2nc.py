from __future__ import division, print_function
import sqlite3
import numpy as np
import xarray as xr
import pandas as pd
from ..core import utils

def virtuosorsk_to_cdf(metadata):
    """
    Main function to load data from RSK file and save to raw .CDF
    """

    ds = virtuosorsk_to_xr(metadata)

    print("Writing to raw netCDF")
    
    # trim
    #ds = utils.clip_ds(ds)

    ds.to_netcdf(ds.attrs['filename'] + '-raw.cdf')

    print("Done")

    return ds, ds.attrs['filename'] + '-raw.cdf'


def init_connection(rskfile):
    """Initialize an sqlite3 connection and return a cursor"""

    conn = sqlite3.connect(rskfile)
    return conn.cursor()

def virtuosorsk_to_xr(metadata):
    """
    Load data from RSK file and generate an xarray Dataset
    
    This version modeled from Dan Nowacki's rsk_to_xr in stglib
    """

    rskfile = metadata['basefile'] + '.rsk'

    ds = xr.Dataset()

    ds = utils.write_metadata(ds, metadata)
    #ds = write_metadata(ds, metadata)

    print(('Loading from sqlite file %s; '
           'this may take a while for large datasets') % rskfile)

    conn = init_connection(rskfile)

    ########  table: deployments  ########
    #['deploymentID', 'instrumentID', 'comment', 'loggerStatus', 'loggerTimeDrift', 'timeOfDownload', 'name', 'sampleSize']
    #(1, 1, '4th attempt\r\nNow with Ruskin v2.3.1', 'stopped', -582653584941, 946684866123, '\\054276_20180618_1614.rsk', 156864)

    # right now we are assuming a database with only on deployment.
    # if there are more than one deployment, this code will work only on 
    # the first deployment for now
    
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
    
    #conn.execute("SELECT tstamp, channel01 FROM burstdata")
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


    ''' # for burst data
    # get indices that end at the end of the final burst
    # datlength = a['unixtime'].shape[0] - a['unixtime'].shape[0] % samplingcount

    # reshape
    for k in a:
        a[k] = a[k][:datlength].reshape(
            (int(datlength/samplingcount), samplingcount)
        )
    times = pd.to_datetime(a['unixtime'][:, 0], unit='ms')
    samples = np.arange(samplingcount)
    '''
    
    # for continuous data
    times = pd.to_datetime(a['unixtime'][:], unit='ms')

    #times = pd.to_datetime(a['unixtime'][:], unit='ms')
    depth = np.array(ds.attrs['WATER_DEPTH']-ds.attrs['initial_instrument_height'])
    lat = ds.attrs['latitude']
    lon = ds.attrs['longitude']
    
    print('we are up to writing Turb')

    ''' # the version for a burst data set, as per dWave
    ds[epic_varname] = xr.DataArray(
        a[varname][:],
        coords=[times, samples],
        dims=('time','sample'),
        name=epic_name,
        attrs={
            'long_name': str(conn.execute(
                "select LongName from channels").fetchall()[0][0]),
            'units': str(conn.execute(
                "select units from channels").fetchall()[0][0]),
            'initial_instrument_height': ds.attrs['initial_instrument_height'],
            'serial_number': ds.attrs['serial_number']},
        encoding={'_FillValue': 1e35})
    '''
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

    # times here are datetime objects <class 'pandas.core.indexes.datetimes.DatetimeIndex'>
    ds['time'] = xr.DataArray(times, dims=('time'), name='time')
    
    # for burst data
    #ds['sample'] = xr.DataArray(samples, dims=('sample'), name = 'sample')
   
    ds['lat'] = xr.DataArray(
        [lat],
        dims=('lat'),
        name='lat',
        attrs={'units': 'degree_north',
               'long_name': 'Latitude',
               'epic_code': 500})

    ds['lon'] = xr.DataArray(
        [lon],
        dims=('lon'),
        name='lon',
        attrs={'units': 'degree_east',
               'long_name': 'Longitude',
               'epic_code': 502})

    ds['depth'] = xr.DataArray(
        [depth],
        dims=('depth'),
        name='depth',
        attrs={'units': 'm',
               'type': 'EVEN',
               'epic_code': 3})

    # need to add  time attrs after DataArrays have been combined into Dataset
    ds['time'].attrs.update({'standard_name': 'time', 'axis': 'T'})

    conn.close()    

    return ds

def virtuosocdf_to_nc(cdf_filename,
              writefile=True,
              format='NETCDF3_64BIT'):
    """
    Load raw .cdf file, trim, apply QAQC, and save to .nc
    """

    # Load raw .cdf data
    ds = xr.open_dataset(cdf_filename, autoclose=True)

    # Clip data to in/out water times or via good_ens
    ds = utils.clip_ds(ds)

    # this function names as epic_time and epic_time2
    ds = utils.create_epic_times(ds)
    
    #ds = create_epic_times_as_primary(ds)

    #ds = utils.create_2d_time(ds)

    ds = Virtuosods_add_attrs(ds)

    # add lat/lon coordinates to each variable
    for var in ds.data_vars:
        if 'time' not in var:
            ds = utils.add_lat_lon(ds, var)
            ds = utils.add_depth_scalar(ds, var)

    ds = utils.add_min_max(ds)

    ds = utils.add_start_stop_time(ds)

    ds = utils.ds_coord_no_fillvalue(ds)

    ds = utils.add_epic_history(ds)

    ds = Virtuoso_add_delta_t(ds)
    
    # not sure if the above operations expect CF time, so do this
    # just before writing
    ds = swap_cf_to_epic_time(ds)

    for var in ds.variables:
        if (var not in ds.coords) and ('time' not in var):
            # cast as float32
            ds = utils.set_var_dtype(ds, var)

    if writefile:
        # Write to .nc file
        print("Writing cleaned/trimmed data to .nc file")
        nc_filename = ds.attrs['filename'] + '-cal.nc'

        ds.to_netcdf(nc_filename, format=format)
        #ds.to_netcdf(nc_filename)

        # Rename time variables for EPIC compliance, keeping a time_cf
        # coorindate.
        if ds.attrs['samples_per_burst'] > 1:
            utils.rename_time_2d(nc_filename)

        print('Done writing netCDF file', nc_filename)

    # rename time variables after the fact to conform with EPIC/CMG standards
    # utils.rename_time(nc_filename)

    return ds, nc_filename

def Virtuosods_add_attrs(ds):
    # Update attributes for EPIC and STG compliance
    ds = utils.ds_coord_no_fillvalue(ds)

    if 'epic_time' in ds:
        ds['time'].attrs.update({'standard_name': 'time',
                                 'axis': 'T'})

        ds['epic_time'].attrs.update({'units': 'True Julian Day',
                                      'type': 'EVEN',
                                      'epic_code': 624})
        
        ds['epic_time2'].attrs.update({'units': 'msec since 0:00 GMT',
                                       'type': 'EVEN',
                                       'epic_code': 624})
    else:
        ds['time'].attrs.update({'units': 'True Julian Day',
                                      'type': 'EVEN',
                                      'epic_code': 624})
    
        ds['time2'].attrs.update({'units': 'msec since 0:00 GMT',
                                       'type': 'EVEN',
                                       'epic_code': 624})

    if 'epic_time_2d' in ds:
        ds['epic_time_2d'].attrs = ds['epic_time'].attrs
    if 'epic_time2_2d' in ds:
        ds['epic_time2_2d'].attrs = ds['epic_time2'].attrs

    ds.attrs['COMPOSITE'] = 0

    return ds

def Virtuoso_add_delta_t(ds):

    ds.attrs['DELTA_T'] = int(ds.attrs['sample_interval'])

    return ds

'''
def create_epic_times_as_primary(ds):
    jd = utils.make_jd(ds['time'].to_dataframe().index)

    ds['time'] = xr.DataArray(utils.make_epic_time(jd),
                                   dims='time',
                                   encoding={'_FillValue': None})

    ds['time2'] = xr.DataArray(utils.make_epic_time2(jd),
                                    dims='time',
                                    encoding={'_FillValue': None})
'''

def swap_cf_to_epic_time(ds):

    if 'epic_time' in ds.variables:
        ds['time'].rename('time_cf')
        ds['epic_time'].rename('time')
        ds['epic_time2'].rename('time2')

    return ds


