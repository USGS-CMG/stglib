#!/usr/bin/env python -m

from __future__ import division, print_function
import sys
import netCDF4
import xarray as xr
from ..core import utils

def nc_to_diwasp(nc_filename):

    ds = xr.open_dataset(nc_filename, autoclose=True, decode_times=False)
    ds['time'] = ds['time_cf']
    ds = ds.drop(['time_cf', 'time2'])
    ds = xr.decode_cf(ds, decode_times=True)

    ds = utils.create_epic_time(ds)

    mat = xr.open_dataset(ds.attrs['filename'][:-2] + 'diwasp.nc', autoclose=True)

    for k in ['wp_peak', 'wh_4061', 'wp_4060']:
        ds[k] = xr.DataArray(mat[k], dims='time')

    ds['frequency'] = xr.DataArray(mat['frequency'], dims=('frequency'))

    ds['pspec'] = xr.DataArray(mat['pspec'], dims=('time', 'frequency'))

    ds = create_water_depth(ds)

    ds = ds.drop(['P_1', 'P_1ac', 'sample'])

    ds = trim_max_wp(ds)

    ds = trim_min_wh(ds)

    ds = trim_wp_ratio(ds)

    # Add attrs
    ds = ds_add_attrs(ds)

    ds = ds_add_diwasp_history(ds)

    nc_filename = ds.attrs['filename'] + 's-a.nc'

    ds = utils.rename_time(ds)

    ds.to_netcdf(nc_filename)

    return ds


def create_water_depth(ds):
    """Create water_depth variable"""

    if 'initial_instrument_height' in ds.attrs:
        if 'P_1ac' in ds:
            ds.attrs['nominal_instrument_depth'] = ds['P_1ac'].mean().values
            ds['water_depth'] = ds.attrs['nominal_instrument_depth']
            wdepth = ds.attrs['nominal_instrument_depth'] + ds.attrs['initial_instrument_height']
            ds.attrs['WATER_DEPTH_source'] = 'water depth = MSL from pressure sensor,'\
                                             ' atmospherically corrected'
            ds.attrs['WATER_DEPTH_datum'] = 'MSL'
        elif 'P_1' in VEL:
            ds.attrs['nominal_instrument_depth'] = ds['P_1'].mean().values
            ds['water_depth'] = ds.attrs['nominal_instrument_depth']
            wdepth = ds.attrs['nominal_instrument_depth'] + ds.attrs['initial_instrument_height']
            ds.attrs['WATER_DEPTH_source'] = 'water depth = MSL from pressure sensor'
            ds.attrs['WATER_DEPTH_datum'] = 'MSL'
        else:
            wdepth = ds.attrs['WATER_DEPTH']
            ds.attrs['nominal_instrument_depth'] = ds.attrs['WATER_DEPTH'] - ds.attrs['initial_instrument_height']
            ds['water_depth'] = ds.attrs['nominal_instrument_depth']
        ds.attrs['WATER_DEPTH'] = wdepth # TODO: why is this being redefined here? Seems redundant
    elif 'nominal_instrument_depth' in ds.attrs:
        ds.attrs['initial_instrument_height'] = ds.attrs['WATER_DEPTH'] - ds.attrs['nominal_instrument_depth']
        ds['water_depth'] = ds.attrs['nominal_instrument_depth']

    if 'initial_instrument_height' not in ds.attrs:
        ds.attrs['initial_instrument_height'] = 0 # TODO: do we really want to set to zero?

    return ds

def trim_max_wp(ds):
    """
    QA/QC
    Trim wave data based on maximum wave period as specified in metadata
    """

    if 'maximum_wp' in ds.attrs:
        print('Trimming using maximum period of %f seconds'
            % ds.attrs['maximum_wp'])
        for var in ['wp_peak', 'wp_4060']:
            ds[var] = ds[var].where((ds['wp_peak'] < ds.attrs['maximum_wp']) &
                (ds['wp_4060'] < ds.attrs['maximum_wp']))

        for var in ['wp_peak', 'wp_4060']:
            notetxt = 'Values filled where wp_peak, wp_4060 >= %f' % ds.attrs['maximum_wp'] + '. '

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

    if 'minimum_wh' in ds.attrs:
        print('Trimming using minimum wave height of %f m'
            % ds.attrs['minimum_wh'])
        ds = ds.where(ds['wh_4061'] > ds.attrs['minimum_wh'])

        for var in ['wp_peak', 'wp_4060', 'wh_4061']:
            notetxt = 'Values filled where wh_4061 <= %f' % ds.attrs['minimum_wh'] + '. '

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
            ds[var] = ds[var].where(ds['wp_peak']/ds['wp_4060'] < ds.attrs['wp_ratio'])

        for var in ['wp_peak', 'wp_4060']:
            notetxt = 'Values filled where wp_peak:wp_4060 >= %f' % ds.attrs['wp_ratio'] + '. '

            if 'note' in ds[var].attrs:
                ds[var].attrs['note'] = notetxt + ds[var].attrs['note']
            else:
                ds[var].attrs.update({'note': notetxt})

    return ds


def ds_add_diwasp_history(ds):
    """
    Add history indicating DIWASP has been applied
    """

    histtext = 'Wave statistics computed using DIWASP 1.4. '

    if 'history' in ds.attrs:
        ds.attrs['history'] = histtext + ds.attrs['history']
    else:
        ds.attrs['history'] = histtext

    return ds

def ds_add_attrs(ds):
    """
    Add EPIC and other attributes to variables
    """

    # Update attributes for EPIC and STG compliance
    ds.lat.encoding['_FillValue'] = False
    ds.lon.encoding['_FillValue'] = False
    ds.depth.encoding['_FillValue'] = False
    ds.time.encoding['_FillValue'] = False
    ds.epic_time.encoding['_FillValue'] = False
    ds.epic_time2.encoding['_FillValue'] = False
    ds.frequency.encoding['_FillValue'] = False

    ds['time'].attrs.update({'standard_name': 'time',
        'axis': 'T'})

    ds['epic_time'].attrs.update({'units': 'True Julian Day',
        'type': 'EVEN',
        'epic_code': 624})

    ds['epic_time2'].attrs.update({'units': 'msec since 0:00 GMT',
        'type': 'EVEN',
        'epic_code': 624})

    def add_attributes(var, dsattrs):
        var.attrs.update({'serial_number': dsattrs['serial_number'],
            'initial_instrument_height': dsattrs['initial_instrument_height'],
            # 'nominal_instrument_depth': metadata['nominal_instrument_depth'], # FIXME
            'height_depth_units': 'm',
            'sensor_type': dsattrs['INST_TYPE'],
            '_FillValue': 1e35})

    ds['wp_peak'].attrs.update({'long_name': 'Dominant (peak) wave period',
        'units': 's',
        'epic_code': 4063})

    ds['wp_4060'].attrs.update({'long_name': 'Average wave period',
        'units': 's',
        'epic_code': 4060})

    ds['wh_4061'].attrs.update({'long_name': 'Significant wave height',
        'units': 'm',
        'epic_code': 4061})

    ds['pspec'].attrs.update({'long_name': 'Pressure derived non-directional wave energy spectrum',
        'units': 'm^2/Hz',
        'note': 'Use caution: all spectra are provisional'})

    ds['frequency'].attrs.update({'long_name': 'Frequency',
        'units': 'Hz'})

    for var in ['wp_peak', 'wh_4061', 'wp_4060', 'pspec', 'water_depth']:
        add_attributes(ds[var], ds.attrs)
        ds[var].attrs.update({'minimum': ds[var].min().values,
            'maximum': ds[var].max().values})

    return ds
