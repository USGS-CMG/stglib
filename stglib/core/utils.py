from __future__ import division, print_function
import csv
import os
import sys
import inspect
import platform
import netCDF4
import warnings
import xarray as xr
import numpy as np
import scipy.io as spio
import pandas as pd
import sqlite3
import stglib


def is_cf(ds):
    if ('Conventions' in ds.attrs) and (str(ds.attrs['Conventions']) == 'CF-1.6'):
       return True
    else:
       return False


def clip_ds(ds, wvs=False):
    """
    Clip an xarray Dataset from metadata, either via good_ens or
    Deployment_date and Recovery_date.

    wvs only applies to Aquadopp waves here. It is distinct from waves flag
    because AQD waves can have a different sampling interval than AQD currents
    """

    print('first burst in full file: {}, idx {}'.format(
        ds['time'].min().values, ds['time'].argmin().values))
    print('last burst in full file: {}, idx {}'.format(
        ds['time'].max().values, ds['time'].argmax().values))

    # clip either by ensemble indices or by the deployment and recovery
    # date specified in metadata
    if 'good_ens' in ds.attrs and not wvs:
        # we have good ensemble indices in the metadata
        print('Clipping data using good_ens')

        # so we can deal with multiple good_ens ranges, or just a single range
        if np.ndim(ds.attrs['good_ens']) == 1:
            good_ens = [ds.attrs['good_ens']]
        else:
            good_ens = ds.attrs['good_ens']

        goods = []
        for x in good_ens:
            goods.append(np.arange(x[0], x[1]))
        goods = np.hstack(goods)

        ds = ds.isel(time=goods)

        histtext = 'Data clipped using good_ens values of {}. '.format(
            str(good_ens))

        ds = insert_history(ds, histtext)

    elif 'good_ens_wvs' in ds.attrs and wvs:
        print('Clipping data using good_ens_wvs')
        good_ens = ds.attrs['good_ens_wvs']

        goods = np.arange(good_ens[0], good_ens[1])

        ds = ds.isel(time=goods)

        histtext = 'Data clipped using good_ens_wvs values of {}. '.format(
            str(good_ens))

        ds = insert_history(ds, histtext)

    elif 'good_dates' in ds.attrs:
        # clip by start/end dates that are not Deployment_date
        # and Recovery_date
        print('Clipping data using good_dates')

        where = np.where(
            (ds['time'].values >= np.datetime64(ds.attrs['good_dates'][0])) &
            (ds['time'].values <= np.datetime64(ds.attrs['good_dates'][1])))[0]
        print('good_dates[0] {}, idx {}'.format(
            ds.attrs['good_dates'][0], where.min()))
        print('good_dates[1] {}, idx {}'.format(
            ds.attrs['good_dates'][1], where.max()))
        ds = ds.sel(time=slice(ds.attrs['good_dates'][0],
                               ds.attrs['good_dates'][1]))

        histtext = 'Data clipped using good_dates of {}. '.format(
            ds.attrs['good_dates'])

        ds = insert_history(ds, histtext)

    elif 'Deployment_date' in ds.attrs and 'Recovery_date' in ds.attrs:
        # we clip by the times in/out of water as specified in the metadata
        print('Clipping data using Deployment_date and Recovery_date')
        ds = ds.sel(time=slice(ds.attrs['Deployment_date'],
                               ds.attrs['Recovery_date']))

        histtext = ('Data clipped using Deployment_date of {} and '
                    'Recovery_date of {}. ').format(
                    ds.attrs['Deployment_date'], ds.attrs['Recovery_date'])

        ds = insert_history(ds, histtext)
    else:
        # do nothing
        print('Did not clip data; no values specified in metadata')

    try:
        print('first burst in trimmed file:', ds['time'].min().values)
        print('last burst in trimmed file:', ds['time'].max().values)
    except ValueError:
        raise(ValueError('No valid time values in trimmed dataset. Are you '
                         'sure you sure you specified Deployment and Recovery '
                         'dates correctly?'))

    return ds


def add_min_max(ds):
    """
    Add minimum and maximum values to variables in NC or CDF files
    This function assumes the data are in xarray DataArrays within Datasets
    """

    exclude = list(ds.dims)
    [exclude.append(k) for k in ds.variables if 'time' in k]
    exclude.extend(['TIM', 'TransMatrix'])

    alloweddims = ['time', 'sample', 'depth']

    for k in ds.variables:
        if k not in exclude:
            kwargs = {'dim': tuple(d for d in alloweddims if d in ds[k].dims)}

            ds[k].attrs.update({
                'minimum': ds[k].min(**kwargs).squeeze().values,
                'maximum': ds[k].max(**kwargs).squeeze().values
                })

    return ds


def insert_history(ds, histtext):
    if 'history' in ds.attrs:
        ds.attrs['history'] = histtext + ds.attrs['history']
    else:
        ds.attrs['history'] = histtext

    return ds


def add_history(ds):
    if is_cf(ds):
        histtext = 'Processed to {} using {}. '.format(
            ds.attrs['Conventions'],
            os.path.basename(sys.argv[0]))
    else:
        histtext = 'Processed to EPIC using {}. '.format(
            os.path.basename(sys.argv[0]))

    return insert_history(ds, histtext)


def ds_add_waves_history(ds):
    histtext = ('Wave statistics computed using scipy.signal.welch(), '
                'assigning cutoff following Jones & Monismith (2007), and '
                'applying f^-4 tail past cutoff. ')

    return insert_history(ds, histtext)


def ds_add_diwasp_history(ds):
    """
    Add history indicating DIWASP has been applied
    """

    histtext = 'Wave statistics computed using DIWASP 1.1GD. '

    return insert_history(ds, histtext)


def ds_coord_no_fillvalue(ds):

    for var in ['lat',
                'lon',
                'depth',
                'time',
                'time2',
                'time_cf',
                'time_2d',
                'time_cf_2d',
                'epic_time',
                'epic_time2',
                'epic_time_2d',
                'epic_time2_2d',
                'TransMatrix',
                'direction',
                'sample',
                'frequency']:
        if var in ds:
            ds[var].encoding['_FillValue'] = None

    return ds


def ds_add_attrs(ds):
    """
    Add EPIC and other attributes to variables
    """

    # Update attributes for EPIC and STG compliance
    ds = ds_coord_no_fillvalue(ds)

    ds['time'].attrs.update({
        'standard_name': 'time',
        'axis': 'T'})

    ds['time'].encoding['dtype'] = 'i4'

    if 'epic_time' in ds:
        ds['epic_time'].attrs.update({
            'units': 'True Julian Day',
            'type': 'EVEN',
            'epic_code': 624})

    if 'epic_time2' in ds:
        ds['epic_time2'].attrs.update({
            'units': 'msec since 0:00 GMT',
            'type': 'EVEN',
            'epic_code': 624})

    def add_attributes(var, dsattrs):
        var.attrs.update({
            'serial_number': dsattrs['serial_number'],
            'initial_instrument_height': dsattrs['initial_instrument_height'],
            'nominal_instrument_depth': dsattrs['nominal_instrument_depth'],
            'height_depth_units': 'm',
            'sensor_type': dsattrs['INST_TYPE']})
        var.encoding['_FillValue'] = 1e35

    ds['wp_peak'].attrs.update({
        'long_name': 'Dominant (peak) wave period',
        'units': 's',
        'epic_code': 4063})

    ds['wp_4060'].attrs.update({
        'long_name': 'Average wave period',
        'units': 's',
        'epic_code': 4060})

    ds['wh_4061'].attrs.update({
        'long_name': 'Significant wave height',
        'units': 'm',
        'epic_code': 4061})

    ds['pspec'].attrs.update({
        'long_name': 'Pressure derived non-directional wave energy spectrum',
        'units': 'm^2/Hz',
        'note': 'Use caution: all spectra are provisional'})

    ds['frequency'].attrs.update({
        'long_name': 'Frequency',
        'units': 'Hz'})

    if 'direction' in ds.coords:
        ds['direction'].attrs.update({
            'long_name': 'Direction (from, relative to true north)',
            'units': 'degrees'})

    if 'dspec' in ds.data_vars:
        ds['dspec'].attrs.update({
            'long_name': 'Directional wave energy spectrum',
            'units': 'm^2/Hz/degree',
            'note': 'Use caution: all spectra are provisional'})

    if 'wvdir' in ds.data_vars:
        ds['wvdir'].attrs.update({
            'long_name': ('Direction of peak period '
                          '(from, relative to true north)'),
            'units': 'degrees',
            'note': ('Compass direction from which waves are propagating as '
                     'defined by the direction with the greatest energy at '
                     'the peak period')})

    if 'dwvdir' in ds.data_vars:
        ds['dwvdir'].attrs.update({
            'long_name': ('Dominant wave direction '
                          '(from, relative to true north)'),
            'units': 'degrees',
            'note': ('Compass direction from which waves are propagating as '
                     'defined by the direction band with greatest total '
                     'energy summed over all frequencies')})

    if 'wd_4062' in ds.data_vars:
        ds['wd_4062'].attrs.update({
            'long_name': 'Mean wave direction',
            'units': 'degrees',
            'epic_code': 4062,
            'note': 'Compass direction from which waves are propagating'})

    for var in ['wp_peak', 'wh_4061', 'wp_4060', 'wd_4062',
                'pspec', 'water_depth', 'dspec', 'wvdir', 'dwvdir']:
        if var in ds.variables:
            add_attributes(ds[var], ds.attrs)
            ds[var].attrs.update({
                'minimum': ds[var].min().values,
                'maximum': ds[var].max().values})

    if 'burst' in ds:
        ds['burst'].encoding['_FillValue'] = 1e35

    return ds


def trim_max_wp(ds):
    """
    QA/QC
    Trim wave data based on maximum wave period as specified in metadata
    """

    if 'wp_max' in ds.attrs:
        print('Trimming using maximum period of {} seconds'.format(
            ds.attrs['wp_max']))
        for var in ['wp_peak', 'wp_4060']:
            ds[var] = ds[var].where(
                    (ds['wp_peak'] < ds.attrs['wp_max']) &
                    (ds['wp_4060'] < ds.attrs['wp_max'])
                    )

        for var in ['wp_peak', 'wp_4060']:
            notetxt = 'Values filled where wp_peak, wp_4060 >= {}. '.format(
                ds.attrs['wp_max'])

            if 'note' in ds[var].attrs:
                ds[var].attrs['note'] = notetxt + ds[var].attrs['note']
            else:
                ds[var].attrs.update({'note': notetxt})

    return ds


def trim_min_wh(ds):
    """
    QA/QC
    Trim wave data based on minimum wave height as specified in metadata
    """

    if 'wh_min' in ds.attrs:
        print('Trimming using minimum wave height of {} m'.format(
            ds.attrs['wh_min']))
        for var in ['wp_peak', 'wh_4061', 'wp_4060']:
            ds[var] = ds[var].where(ds['wh_4061'] > ds.attrs['wh_min'])

            notetxt = 'Values filled where wh_4061 <= {}. '.format(
                ds.attrs['wh_min'])

            if 'note' in ds[var].attrs:
                ds[var].attrs['note'] = notetxt + ds[var].attrs['note']
            else:
                ds[var].attrs.update({'note': notetxt})

    return ds


def trim_max_wh(ds):
    """
    QA/QC
    Trim wave data based on maximum wave height as specified in metadata
    """

    if 'wh_max' in ds.attrs:
        print('Trimming using maximum wave height of {} m'.format(
            ds.attrs['wh_max']))
        for var in ['wp_peak', 'wh_4061', 'wp_4060']:
            ds[var] = ds[var].where(ds['wh_4061'] < ds.attrs['wh_max'])

            notetxt = 'Values filled where wh_4061 >= {}. '.format(
                ds.attrs['wh_max'])

            if 'note' in ds[var].attrs:
                ds[var].attrs['note'] = notetxt + ds[var].attrs['note']
            else:
                ds[var].attrs.update({'note': notetxt})

    return ds


def trim_wp_ratio(ds):
    """
    QA/QC
    Trim wave data based on maximum ratio of wp_peak to wp_4060
    """

    if 'wp_ratio' in ds.attrs:
        print('Trimming using maximum ratio of wp_peak to wp_4060 of %f'
              % ds.attrs['wp_ratio'])
        for var in ['wp_peak', 'wp_4060']:
            ds[var] = ds[var].where(
                ds['wp_peak']/ds['wp_4060'] < ds.attrs['wp_ratio'])

        for var in ['wp_peak', 'wp_4060']:
            notetxt = 'Values filled where wp_peak:wp_4060 >= {}. '.format(
                ds.attrs['wp_ratio'])

            if 'note' in ds[var].attrs:
                ds[var].attrs['note'] = notetxt + ds[var].attrs['note']
            else:
                ds[var].attrs.update({'note': notetxt})

    return ds


def write_metadata(ds, metadata):
    """Write out all metadata to CDF file"""

    for k in metadata:
        # recursively write out instmeta
        if k == 'instmeta':
            ds.attrs.update({x: metadata[k][x] for x in metadata[k]})
        else:
            ds.attrs.update({k: metadata[k]})

    f = os.path.basename(inspect.stack()[1][1])

    histtext = ('Processed using {} with stglib {}, xarray {}, NumPy {}, '
                'netCDF4 {}, Python {}. ').format(
                    f, stglib.__version__, xr.__version__, np.__version__,
                    netCDF4.__version__, platform.python_version())

    ds = insert_history(ds, histtext)

    return ds


def set_var_dtype(ds, var, dtype='float32'):

    ds[var].encoding['dtype'] = dtype

    return ds


def rename_time(ds):
    """
    Rename time variables for EPIC compliance, keeping a time_cf coorindate.
    """

    if is_cf(ds):
        pass
    else:
        ds = ds.rename({'time': 'time_cf'})
        ds = ds.rename({'epic_time': 'time'})
        ds = ds.rename({'epic_time2': 'time2'})
        ds = ds.set_coords(['time', 'time2'])
        ds = ds.swap_dims({'time_cf': 'time'})
        # output int32 time_cf for THREDDS compatibility
        ds['time_cf'].encoding['dtype'] = 'i4'

    return ds


def rename_time_2d(nc_filename, ds):
    if is_cf(ds):
        print('not renaming 2D time because CF==1.6')
        pass
    else:
        # Need to do this in two steps after renaming the variable.
        # Not sure why, but it works this way.
        with netCDF4.Dataset(nc_filename, 'r+') as nc:
            nc.renameVariable('time', 'time_cf')
            # nc.renameVariable('time_2d', 'time_cf_2d')
            timebak = nc['epic_time_2d'][:]
            nc.renameVariable('epic_time_2d', 'time')
            nc.renameVariable('epic_time2_2d', 'time2')

        with netCDF4.Dataset(nc_filename, 'r+') as nc:
            nc['time'][:] = timebak


def open_time_2d_dataset(filename):
    # need to drop 'time' variable because of xarray limitations related
    # to coordinates and variables with the same name, otherwise it raises a
    # MissingDimensionsError
    # Check if CF or not, and return the correct dataset
    with xr.open_dataset(filename,
                         decode_times=False,
                         drop_variables='time') as ds:
        if is_cf(ds):
            iscf = True
        else:
            iscf = False

    if iscf:
        return xr.open_dataset(filename)
    else:
        return xr.open_dataset(filename,
                               decode_times=False,
                               drop_variables='time')


def epic_to_cf_time(ds):
    if 'time_cf' in ds:
        ds['time'] = ds['time_cf']
    else:
        ds['time'] = epic_to_datetime(ds['time'].values, ds['time2'].values)

    for v in ['time_cf', 'time2', 'epic_time', 'epic_time2']:
        if v in ds:
            ds = ds.drop(v)
    return xr.decode_cf(ds, decode_times=True)


def epic_to_datetime(time, time2):
    thedate = pd.to_datetime(time - 0.5, origin='julian', unit='D')
    thetime = pd.to_timedelta(time2, unit='ms')
    return thedate + thetime


def make_jd(time):
    return time.to_julian_date().values + 0.5


def make_epic_time(jd):
    epic_time = np.floor(jd)
    # make sure they are all integers, and then cast as such
    if np.all(np.mod(epic_time, 1) == 0):
        epic_time = epic_time.astype(np.int32)
    else:
        warnings.warn('Not all EPIC time values are integers; '
                      'this will cause problems with time and time2')

    return epic_time


def make_epic_time2(jd):
    # TODO: Hopefully this is correct... roundoff errors on big numbers...
    return np.round((jd - np.floor(jd))*86400000).astype(np.int32)


def create_epic_times(ds, waves=False):
    jd = make_jd(ds['time'].to_dataframe().index)

    ds['epic_time'] = xr.DataArray(make_epic_time(jd), dims='time')
    ds['epic_time'].encoding['_FillValue'] = None

    ds['epic_time2'] = xr.DataArray(make_epic_time2(jd), dims='time')
    ds['epic_time2'].encoding['_FillValue'] = None

    return ds


def create_2d_time(ds):
    print('Creating 2D time variable')
    # time increment in milliseconds
    td = (ds.attrs['sample_interval'] *
          np.arange(ds.attrs['samples_per_burst']) * 1000)

    # time_2d is a CF representation of a 2d time
    ds['time_2d'] = xr.DataArray(np.expand_dims(ds['time'], 1) +
                                 [np.timedelta64(int(x), 'ms') for x in td],
                                 dims=('time', 'sample'))

    raveljd = make_jd(pd.DatetimeIndex(np.ravel(ds['time_2d'])))
    jd_2d = np.reshape(raveljd, ds['time_2d'].shape)

    ds['epic_time_2d'] = xr.DataArray(make_epic_time(jd_2d),
                                      dims=('time', 'sample'))
    ds['epic_time_2d'].encoding['_FillValue'] = None

    ds['epic_time2_2d'] = xr.DataArray(make_epic_time2(jd_2d),
                                       dims=('time', 'sample'))
    ds['epic_time2_2d'].encoding['_FillValue'] = None

    ds = ds.drop('time_2d')  # don't need it anymore

    return ds


def add_start_stop_time(ds):
    """Add start_time and stop_time attrs"""

    ds.attrs.update({'start_time': ds['time'][0].values.astype(str),
                     'stop_time': ds['time'][-1].values.astype(str)})

    return ds


def add_lat_lon(ds, var):
    """Add lat and lon dimensions"""

    ds[var] = xr.concat([ds[var]], dim=ds['lon'])
    ds[var] = xr.concat([ds[var]], dim=ds['lat'])

    # Reorder so lat, lon are at the end.
    dims = [d for d in ds[var].dims if (d != 'lon') and (d != 'lat')]
    dims.extend(['lat', 'lon'])
    dims = tuple(dims)

    ds[var] = ds[var].transpose(*dims)

    return ds


def ds_add_lat_lon(ds):
    ds['lat'] = xr.DataArray(
        [ds.attrs['latitude']],
        dims=('lat'),
        name='lat',
        attrs={'units': 'degree_north',
               'long_name': 'Latitude',
               'standard_name': 'latitude',
               'epic_code': 500})

    ds['lon'] = xr.DataArray(
        [ds.attrs['longitude']],
        dims=('lon'),
        name='lon',
        attrs={'units': 'degree_east',
               'long_name': 'Longitude',
               'standard_name': 'longitude',
               'epic_code': 502})

    return ds


def shift_time(ds, timeshift):
    """Shift time to middle of burst"""

    if timeshift != 0:
        # shift times to center of ensemble
        ds['time'] = ds['time'] + np.timedelta64(int(timeshift), 's')
        print('Time shifted by %.f s' % int(timeshift))
        if not timeshift.is_integer():
            warnings.warn(
                'time offset of %.3f s was adjusted to %.f s for shifting time' %
                (timeshift, int(timeshift)))

    if 'ClockError' in ds.attrs:
        if ds.attrs['ClockError'] != 0:
            # note negative on ds.attrs['ClockError']
            ds['time'] = (ds['time'] +
                          np.timedelta64(-ds.attrs['ClockError'], 's'))

        histtext = 'Time shifted by {:d} s from ClockError. '.format(
            -ds.attrs['ClockError'])

        insert_history(ds, histtext)

    return ds


def create_water_depth_var(ds):

    press = None

    if 'Pressure_ac' in ds:
        press = 'Pressure_ac'
    elif 'P_1ac' in ds:
        press = 'P_1ac'
    elif 'Pressure' in ds:
        press = 'Pressure'
    elif 'P_1' in ds:
        press = 'P_1'

    if 'sample' in ds.dims:
        ds['water_depth'] = xr.DataArray(
            ds[press].squeeze().mean(dim='sample') +
            ds.attrs['initial_instrument_height'])
    else:
        ds['water_depth'] = xr.DataArray(
            ds[press].squeeze() + ds.attrs['initial_instrument_height'])

    ds['water_depth'].attrs['long_name'] = 'Total water depth'
    ds['water_depth'].attrs['units'] = 'm'

    return ds


def create_water_depth(ds):
    """Create WATER_DEPTH attribute"""

    press = None

    if 'Pressure_ac' in ds:
        press = 'Pressure_ac'
    elif 'P_1ac' in ds:
        press = 'P_1ac'
    elif 'Pressure' in ds:
        press = 'Pressure'
    elif 'P_1' in ds:
        press = 'P_1'

    if 'sample' in ds.dims:
        dims = ('time', 'sample')
    else:
        dims = 'time'

    if 'initial_instrument_height' in ds.attrs:
        if press:
            ds.attrs['nominal_instrument_depth'] = (
                ds[press].squeeze().mean(dim=dims).values)
            # ds['water_depth'] = ds.attrs['nominal_instrument_depth']
            wdepth = (
                ds.attrs['nominal_instrument_depth'] +
                ds.attrs['initial_instrument_height'])
            if 'ac' in press:
                ds.attrs['WATER_DEPTH_source'] = ('water depth = MSL from '
                                                  'pressure sensor, '
                                                  'atmospherically corrected')
            else:
                ds.attrs['WATER_DEPTH_source'] = ('water depth = MSL from '
                                                  'pressure sensor')
            ds.attrs['WATER_DEPTH_datum'] = 'MSL'
        else:
            wdepth = ds.attrs['WATER_DEPTH']
            ds.attrs['nominal_instrument_depth'] = (
                ds.attrs['WATER_DEPTH'] -
                ds.attrs['initial_instrument_height'])
        # ds['Depth'] = ds.attrs['nominal_instrument_depth']
        # TODO: why is this being redefined here? Seems redundant
        ds.attrs['WATER_DEPTH'] = wdepth

    elif 'nominal_instrument_depth' in ds.attrs:
        ds.attrs['initial_instrument_height'] = (
            ds.attrs['WATER_DEPTH'] -
            ds.attrs['nominal_instrument_depth'])
        # ds['water_depth'] = ds.attrs['nominal_instrument_depth']

    if 'initial_instrument_height' not in ds.attrs:
        # TODO: do we really want to set to zero?
        ds.attrs['initial_instrument_height'] = 0

    return ds


def create_nominal_instrument_depth(ds):
    if 'nominal_instrument_depth' not in ds.attrs:
        ds.attrs['nominal_instrument_depth'] = (
            ds.attrs['WATER_DEPTH'] - ds.attrs['initial_instrument_height'])

    return ds


def no_p_create_depth(ds):
    # no_p = no pressure sensor. also use for exo
    if 'NAVD88_ref' in ds.attrs:
        ds['depth'] = xr.DataArray(
            [-ds.attrs['NAVD88_ref'] - ds.attrs['initial_instrument_height']],
            dims='depth')
        ds['depth'].attrs['VERT_DATUM'] = 'NAVD88'
        ds['depth'].attrs['NOTE'] = ('Computed as platform depth '
                                     '[m NAVD88] minus '
                                     'initial_instrument_height')
    else:
        ds['depth'] = xr.DataArray(
            [ds.attrs['WATER_DEPTH'] - ds.attrs['initial_instrument_height']],
            dims='depth')
        ds['depth'].attrs['NOTE'] = ('Computed as WATER_DEPTH minus '
                                     'initial_instrument_height')

    ds['depth'].attrs['positive'] = 'down'
    ds['depth'].attrs['axis'] = 'z'
    ds['depth'].attrs['units'] = 'm'
    ds['depth'].attrs['epic_code'] = 3
    ds['depth'].encoding['_FillValue'] = None

    return ds


def no_p_add_depth(ds, var):
    # no_p = no pressure sensor. also use for exo
    ds[var] = xr.concat([ds[var]], dim=ds['depth'])

    # Reorder so lat, lon are at the end.
    dims = [d for d in ds[var].dims if (d != 'depth')]
    dims.extend(['depth'])
    dims = tuple(dims)

    ds[var] = ds[var].transpose(*dims)

    return ds


def insert_note(ds, var, notetxt):
    if 'note' in ds[var].attrs:
        ds[var].attrs['note'] = notetxt + ds[var].attrs['note']
    else:
        ds[var].attrs.update({'note': notetxt})

    return ds


def add_delta_t(ds):
    deltat = np.asscalar(
        (ds['time'][1] - ds['time'][0]) / np.timedelta64(1, 's'))
    if not deltat.is_integer():
        warnings.warn('DELTA_T is not an integer; casting as int in attrs')

    ds.attrs['DELTA_T'] = int(deltat)

    return ds


def check_valid_metadata(metadata):
    for k in ['initial_instrument_height',
              'orientation']:
        if k not in metadata:
            raise KeyError(
                k + ' must be defined, most likely in config.yaml')


def read_samplingrates_burst(ds, conn):
    '''
    Reads in sample information from RBR instrument in burst mode
    '''
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

    # Get sampling interval
    try:
        samplingperiod = conn.execute(
            "select samplingperiod from schedules").fetchall()[0][0]
    except sqlite3.OperationalError:
        samplingperiod = conn.execute(
            "select samplingperiod from wave").fetchall()[0][0]
    ds.attrs['sample_interval'] = samplingperiod / 1000

    # Get Repetition interval of sampling
    try:
        repetitionperiod = conn.execute(
            "select repetitionperiod from schedules").fetchall()[0][0]
    except sqlite3.OperationalError:
        repetitionperiod = conn.execute(
            "select repetitionperiod from wave").fetchall()[0][0]

    # Convert to seconds
    ds.attrs['burst_interval'] = repetitionperiod / 1000

    # Length of bursts in data points
    ds.attrs['burst_length'] = ds.attrs['samples_per_burst'] * \
        ds.attrs['sample_interval']

    return ds


def read_samplingrates_continuous(ds, conn):
    '''
    Reads in sample information from RBR instrument in continuous mode
    '''
    try:
        samplingperiod = conn.execute(
            "select samplingperiod from schedules").fetchall()[0][0]
    except sqlite3.OperationalError:
        samplingperiod = conn.execute(
            "select samplingperiod from continuous").fetchall()[0][0]

    samplingperiod = samplingperiod / 1000  # convert from ms to sec
    samplingrate = 1 / samplingperiod  # convert to rate [Hz]

    # Set sampling period, [sec]
    ds.attrs['sample_interval'] = samplingperiod

    # Set samples per burst
    samplingcount = ds.attrs['wave_interval']*samplingrate
    ds.attrs['samples_per_burst'] = round(samplingcount)

    # Set burst interval, [sec], USER DEFINED in instrument attr
    ds.attrs['burst_interval'] = ds.attrs['wave_interval']

    # Set sample interval
    ds.attrs['burst_length'] = ds.attrs['samples_per_burst'] * \
        ds.attrs['sample_interval']

    return ds


def read_globalatts(fname):
    """
    Read global attributes file (glob_attxxxx.txt) and create metadata
    structure
    """

    metadata = {}

    with open(fname, 'r') as csvfile:
        a = csv.reader(csvfile, delimiter=';')

        for row in a:
            if row[0].strip() == 'MOORING':
                metadata[row[0].strip()] = str(row[1].strip())
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

    from: `StackOverflow <https://stackoverflow.com/q/7008608>`_
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
