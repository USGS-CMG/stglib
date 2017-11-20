#!/usr/bin/env python

from __future__ import division, print_function

import sqlite3
import inspect
import platform
import os
import numpy as np
import xarray as xr
import pandas as pd
import netCDF4


def rsk_to_cdf(metadata):
    """
    Main function to load data from RSK file and save to raw .CDF
    """

    RAW, metadata = rsk_to_xr(metadata)

    print("Writing to raw netCDF")
    xr_to_cdf(RAW, metadata)

    print("Done")

    return RAW, metadata

def init_connection(rskfile):
    """Initialize an sqlite3 connection and return a cursor"""

    conn = sqlite3.connect(rskfile)
    return conn.cursor()

def rsk_to_xr(metadata):
    """
    Load data from RSK file and generate an xarray Dataset
    """

    rskfile = metadata['basefile'] + '.rsk'

    print('Loading from sqlite file %s; this may take a while for large datasets' % rskfile)

    conn = init_connection(rskfile)

    conn.execute("SELECT tstamp, channel01 FROM burstdata")
    data = conn.fetchall()
    print("Done fetching data")
    d = np.asarray(data)
    # Get samples per burst
    samplingcount = conn.execute("select samplingcount from schedules").fetchall()[0][0]

    metadata['samples_per_burst'] = samplingcount
    samplingperiod = conn.execute("select samplingperiod from schedules").fetchall()[0][0]
    metadata['sample_interval'] = samplingperiod / 1000
    repetitionperiod = conn.execute("select repetitionperiod from schedules").fetchall()[0][0]
    metadata['burst_interval'] = repetitionperiod / 1000
    metadata['burst_length'] = metadata['samples_per_burst'] * metadata['sample_interval']
    metadata['serial_number'] = conn.execute("select serialID from instruments").fetchall()[0][0]
    metadata['INST_TYPE'] = 'RBR Virtuoso d|wave'

    a = {}
    a['unixtime'] = d[:,0].copy()
    a['pres'] = d[:,1].copy()
    # sort by time (not sorted for some reason)
    sort = np.argsort(a['unixtime'])
    a['unixtime'] = a['unixtime'][sort]
    a['pres'] = a['pres'][sort]

    # get indices that end at the end of the final burst
    datlength = a['unixtime'].shape[0] - a['unixtime'].shape[0] % samplingcount

    # reshape
    for k in a:
        a[k] = a[k][:datlength].reshape((int(datlength/samplingcount), samplingcount))

    times = pd.to_datetime(a['unixtime'][:,0], unit='ms')
    samples = np.arange(samplingcount)

    dwave = {}

    dwave['P_1'] = xr.DataArray(a['pres'], coords=[times, samples],
        dims=('time', 'sample'), name='Pressure',
        attrs={'long_name': 'Pressure',
               '_FillValue': 1e35,
               'units': 'dbar',
               'epic_code': 1,
               'height_depth_units': 'm',
               'initial_instrument_height': metadata['initial_instrument_height'],
               'serial_number': metadata['serial_number']})

    dwave['time'] = xr.DataArray(times, dims=('time'), name='time')

    dwave['sample'] = xr.DataArray(samples, dims=('sample'), name='sample')

    dwave['lat'] = xr.DataArray([metadata['latitude']], dims=('lat'), name='lat',
        attrs={'units': 'degree_north',
               'long_name': 'Latitude',
               'epic_code': 500})

    dwave['lon'] = xr.DataArray([metadata['longitude']], dims=('lon'), name='lon',
        attrs={'units': 'degree_east',
               'long_name': 'Longitude',
               'epic_code': 502})

    dwave['depth'] = xr.DataArray([metadata['WATER_DEPTH']], dims=('depth'), name='depth',
        attrs={'units': 'm',
               'long_name': 'mean water depth',
               'axis': 'z', # TODO: are these attrs necessary/appropriate?
               'positive': 'down', # TODO: are these attrs necessary/appropriate?
               'epic_code': 3})

    # Create Dataset from dictionary of DataArrays
    RAW = xr.Dataset(dwave)

    # need to add the time attrs after DataArrays have been combined into Dataset
    RAW['time'].attrs.update({'standard_name': 'time', 'axis': 'T'})

    RAW = write_metadata(RAW, metadata)

    return RAW, metadata

def xr_to_cdf(RAW, metadata):
    """Write raw xarray Dataset to .cdf"""

    cdf_filename = metadata['filename'] + '-raw.cdf'

    RAW.to_netcdf(cdf_filename)

def write_metadata(ds, metadata):
    """Write metadata to Dataset"""

    for k in metadata:
        ds.attrs.update({k: metadata[k]})

    f = os.path.basename(inspect.stack()[1][1])

    ds.attrs.update({'history': 'Processed using ' + f + ' with Python ' +
        platform.python_version() + ', xarray ' + xr.__version__ + ', NumPy ' +
        np.__version__ + ', netCDF4 ' + netCDF4.__version__})

    return ds

        # # TODO: add the following??
        # # {'positive','down';
        # #                'long_name', 'Depth';
        # #                'axis','z';
        # #                'units', 'm';
        # #                'epic_code', 3};
        #
        # Pressid = rg.createVariable('Pressure', 'f', ('time','sample',), fill_value=False, zlib=True)
        # Pressid.units = 'dbar'
        # Pressid.long_name = 'Pressure (dbar)'
        # Pressid.generic_name = 'press'
        # Pressid.note = 'raw pressure from instrument, not corrected for changes in atmospheric pressure'
        # Pressid.epic_code = 1
        # Pressid.height_depth_units = 'm'

def main():
    import sys
    sys.path.insert(0, '/Users/dnowacki/Documents/rsklib')
    import rsklib
    import argparse
    import yaml

    parser = argparse.ArgumentParser(description='Convert raw RBR d|wave files (.rsk) to raw .cdf format. Run this script from the directory containing d|wave files')
    parser.add_argument('gatts', help='path to global attributes file (gatts formatted)')
    parser.add_argument('config', help='path to ancillary config file (YAML formatted)')

    args = parser.parse_args()

    # initialize metadata from the globalatts file
    metadata = rsklib.read_globalatts(args.gatts)

    # Add additional metadata from metadata config file
    config = yaml.safe_load(open(args.config))

    for k in config:
        metadata[k] = config[k]

    RAW, metadata = rsklib.rsk_to_cdf(metadata)

    return RAW, metadata

if __name__ == '__main__':
    main()
