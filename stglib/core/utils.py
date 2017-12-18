from __future__ import division, print_function
import csv
import os
import sys
import inspect
import platform
import netCDF4
import xarray as xr
import numpy as np
import scipy.io as spio

def clip_ds(ds):
    """
    Clip an xarray Dataset from metadata, either via good_ens or
    Deployment_date and Recovery_date
    """

    print('first burst in full file:', ds['time'].min().values)
    print('last burst in full file:', ds['time'].max().values)

    # clip either by ensemble indices or by the deployment and recovery date specified in metadata
    if 'good_ens' in ds.attrs:
        # we have good ensemble indices in the metadata
        print('Clipping data using good_ens')

        ds = ds.isel(time=slice(ds.attrs['good_ens'][0], ds.attrs['good_ens'][1]))

        histtext = 'Data clipped using good_ens values of ' + ds.attrs['good_ens'][0] + ', ' + ds.attrs['good_ens'][1] + '. '
        if 'history' in ds.attrs:
            ds.attrs['history'] = histtext + ds.attrs['history']
        else:
            ds.attrs['history'] = histtext

    elif 'Deployment_date' in ds.attrs and 'Recovery_date' in ds.attrs:
        # we clip by the times in/out of water as specified in the metadata
        print('Clipping data using Deployment_date and Recovery_date')

        ds = ds.sel(time=slice(ds.attrs['Deployment_date'], ds.attrs['Recovery_date']))

        histtext = 'Data clipped using Deployment_date and Recovery_date of ' + ds.attrs['Deployment_date'] + ', ' + ds.attrs['Recovery_date'] + '. '
        if 'history' in ds.attrs:
            ds.attrs['history'] = histtext + ds.attrs['history']
        else:
            ds.attrs['history'] = histtext
    else:
        # do nothing
        print('Did not clip data; no values specified in metadata')

    print('first burst in trimmed file:', ds['time'].min().values)
    print('last burst in trimmed file:', ds['time'].max().values)

    return ds

def add_min_max(ds):
    """
    Add minimum and maximum values to variables in NC or CDF files
    This function assumes the data are in xarray DataArrays within Datasets
    """

    exclude = list(ds.dims)
    exclude.extend(('epic_time', 'epic_time2', 'time', 'time2', 'TIM'))

    for k in ds.variables:
        if k not in exclude:
            ds[k].attrs.update({'minimum': ds[k].min().values,
                                'maximum': ds[k].max().values})

    return ds

def add_epic_history(ds):

    ds.attrs['history'] = 'Processed to EPIC using ' + \
                          os.path.basename(sys.argv[0]) + \
                          '. ' + ds.attrs['history']

    return ds


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

def rename_time(ds):
    """
    Rename time variables for EPIC compliance, keeping a time_cf coorindate.
    """

    # nc = netCDF4.Dataset(nc_filename, 'r+')
    # timebak = nc['epic_time'][:]
    # nc.renameVariable('time', 'time_cf')
    # nc.renameVariable('epic_time', 'time')
    # nc.renameVariable('epic_time2', 'time2')
    # nc.close()
    #
    # # need to do this in two steps after renaming the variable
    # # not sure why, but it works this way
    # nc = netCDF4.Dataset(nc_filename, 'r+')
    # nc['time'][:] = timebak
    # nc.close()

    ds.rename({'time': 'time_cf'}, inplace=True)
    ds.rename({'epic_time': 'time'}, inplace=True)
    ds.rename({'epic_time2': 'time2'}, inplace=True)
    ds.set_coords(['time', 'time2'], inplace=True)
    ds.swap_dims({'time_cf': 'time'}, inplace=True)

    return ds

def create_epic_time(ds):

    # create Julian date
    ds['jd'] = ds['time'].to_dataframe().index.to_julian_date() + 0.5

    ds['epic_time'] = np.floor(ds['jd'])
    if np.all(np.mod(ds['epic_time'], 1) == 0): # make sure they are all integers, and then cast as such
        ds['epic_time'] = ds['epic_time'].astype(np.int32)
    else:
        warnings.warn('not all EPIC time values are integers; '\
                      'this will cause problems with time and time2')

    # TODO: Hopefully this is correct... roundoff errors on big numbers...
    ds['epic_time2'] = np.round((ds['jd'] - np.floor(ds['jd']))*86400000).astype(np.int32)

    return ds

def add_start_stop_time(ds):
    """Add start_time and stop_time attrs"""

    ds.attrs.update({'start_time': ds['time'][0].values.astype(str),
                     'stop_time': ds['time'][-1].values.astype(str)})

    return ds

def create_water_depth(ds):
    """Create water_depth variable"""

    if 'initial_instrument_height' in ds.attrs:
        if 'Pressure_ac' in ds:
            ds.attrs['nominal_instrument_depth'] = np.nanmean(ds['Pressure_ac'])
            ds['Depth'] = ds.attrs['nominal_instrument_depth']
            wdepth = ds.attrs['nominal_instrument_depth'] + ds.attrs['initial_instrument_height']
            ds.attrs['WATER_DEPTH_source'] = 'water depth = MSL from pressure sensor, atmospherically corrected'
            ds.attrs['WATER_DEPTH_datum'] = 'MSL'
        elif 'Pressure' in ds:
            ds.attrs['nominal_instrument_depth'] = np.nanmean(ds['Pressure'])
            ds['Depth'] = ds.attrs['nominal_instrument_depth']
            wdepth = ds.attrs['nominal_instrument_depth'] + ds.attrs['initial_instrument_height']
            ds.attrs['WATER_DEPTH_source'] = 'water depth = MSL from pressure sensor'
            ds.attrs['WATER_DEPTH_datum'] = 'MSL'
        else:
            wdepth = ds.attrs['WATER_DEPTH']
            ds.attrs['nominal_instrument_depth'] = ds.attrs['WATER_DEPTH'] - ds.attrs['initial_instrument_height']
            ds['Depth'] = ds.attrs['nominal_instrument_depth']
        ds.attrs['WATER_DEPTH'] = wdepth # TODO: why is this being redefined here? Seems redundant
    elif 'nominal_instrument_depth' in ds.attrs:
        ds.attrs['initial_instrument_height'] = ds.attrs['WATER_DEPTH'] - ds.attrs['nominal_instrument_depth']
        ds['Depth'] = ds.attrs['nominal_instrument_depth']

    if 'initial_instrument_height' not in ds.attrs:
        ds.attrs['initial_instrument_height'] = 0 # TODO: do we really want to set to zero?

    return ds

def read_globalatts(fname):
    """
    Read global attributes file (glob_attxxxx.txt) and create metadata structure
    """

    metadata = {}

    with open(fname, 'r') as csvfile:
        a = csv.reader(csvfile, delimiter=';')

        for row in a:
            if row[0] == 'MOORING':
                metadata[row[0].strip()] = row[1].strip()
            else:
                metadata[row[0].strip()] = str2num(row[1].strip())

        return metadata

def str2num(s):
    """
    Convert string to float if possible
    """

    try:
        float(s)
        return float(s)
    except ValueError:
        return s

def loadmat(filename):
    '''
    this function should be called instead of direct spio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects

    from: `StackOverflow <http://stackoverflow.com/questions/7008608/scipy-io-loadmat-nested-structures-i-e-dictionaries>`_
    '''
    data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)

def _check_keys(dic):
    '''
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries
    '''
    for key in dic:
        if isinstance(dic[key], spio.matlab.mio5_params.mat_struct):
            dic[key] = _todict(dic[key])
    return dic

def _todict(matobj):
    '''
    A recursive function which constructs from matobjects nested dictionaries
    '''
    dic = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, spio.matlab.mio5_params.mat_struct):
            dic[strg] = _todict(elem)
        else:
            dic[strg] = elem
    return dic
