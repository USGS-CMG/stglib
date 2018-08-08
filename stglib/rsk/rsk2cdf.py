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

    rskfile = metadata['basefile'] + '.rsk'

    ds = xr.Dataset()

    ds = utils.write_metadata(ds, metadata)

    print(('Loading from sqlite file %s; '
           'this may take a while for large datasets') % rskfile)

    conn = init_connection(rskfile)

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

    conn.close()

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
