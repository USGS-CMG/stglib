#!/usr/bin/env python

from __future__ import division, print_function
import sys
sys.path.insert(0, '/Users/dnowacki/Documents/rsklib')
import rsklib
import netCDF4
import xarray as xr


def nc_to_diwasp(metadata):

    ds = xr.open_dataset(metadata['filename'] + 'b-cal.nc', autoclose=True, decode_times=False)
    ds['time'] = ds['time_cf']
    ds = ds.drop(['time_cf', 'time2'])
    ds = xr.decode_cf(ds, decode_times=True)

    ds = rsklib.rskcdf2nc.create_epic_time(ds)

    mat = xr.open_dataset(metadata['filename'][:-2] + 'diwasp.nc', autoclose=True)

    for k in ['wp_peak', 'wh_4061', 'wp_4060']:
        ds[k] = xr.DataArray(mat[k], dims='time')

    ds['frequency'] = xr.DataArray(mat['frequency'], dims=('frequency'))

    ds['pspec'] = xr.DataArray(mat['pspec'], dims=('time', 'frequency'))

    ds, metadata = create_water_depth(ds, metadata)

    ds = ds.drop(['P_1', 'P_1ac', 'sample'])

    ds = trim_max_wp(ds, metadata)

    ds = trim_min_wh(ds, metadata)

    ds = trim_wp_ratio(ds, metadata)

    # Add attrs
    ds = ds_add_attrs(ds, metadata)

    ds = rsklib.write_metadata(ds, metadata)

    write_nc(ds, metadata)

    return ds


def write_nc(ds, metadata):
    """Write cleaned and trimmed Dataset to .nc file"""

    nc_filename = metadata['filename'] + 's-a.nc'

    ds.to_netcdf(nc_filename, unlimited_dims='time', engine='netcdf4')

    # rename time variables after the fact to conform with EPIC/CMG standards
    rename_time(nc_filename)


def create_water_depth(VEL, metadata):
    """Create water_depth variable"""

    if 'initial_instrument_height' in metadata:
        if 'P_1ac' in VEL:
            metadata['nominal_instrument_depth'] = VEL['P_1ac'].mean().values
            VEL['water_depth'] = metadata['nominal_instrument_depth']
            wdepth = metadata['nominal_instrument_depth'] + metadata['initial_instrument_height']
            metadata['WATER_DEPTH_source'] = 'water depth = MSL from pressure sensor,'\
                                             ' atmospherically corrected'
            metadata['WATER_DEPTH_datum'] = 'MSL'
        elif 'P_1' in VEL:
            metadata['nominal_instrument_depth'] = VEL['P_1'].mean().values
            VEL['water_depth'] = metadata['nominal_instrument_depth']
            wdepth = metadata['nominal_instrument_depth'] + metadata['initial_instrument_height']
            metadata['WATER_DEPTH_source'] = 'water depth = MSL from pressure sensor'
            metadata['WATER_DEPTH_datum'] = 'MSL'
        else:
            wdepth = metadata['WATER_DEPTH']
            metadata['nominal_instrument_depth'] = metadata['WATER_DEPTH'] - metadata['initial_instrument_height']
            VEL['water_depth'] = metadata['nominal_instrument_depth']
        metadata['WATER_DEPTH'] = wdepth # TODO: why is this being redefined here? Seems redundant
    elif 'nominal_instrument_depth' in metadata:
        metadata['initial_instrument_height'] = metadata['WATER_DEPTH'] - metadata['nominal_instrument_depth']
        VEL['water_depth'] = metadata['nominal_instrument_depth']

    if 'initial_instrument_height' not in metadata:
        metadata['initial_instrument_height'] = 0 # TODO: do we really want to set to zero?

    return VEL, metadata

def trim_max_wp(ds, metadata):
    """
    QA/QC
    Trim wave data based on maximum wave period as specified in metadata
    """

    if 'maximum_wp' in metadata:
        print('Trimming using maximum period of %f seconds'
            % metadata['maximum_wp'])
        for var in ['wp_peak', 'wp_4060']:
            ds[var] = ds[var].where((ds['wp_peak'] < metadata['maximum_wp']) &
                (ds['wp_4060'] < metadata['maximum_wp']))

        for var in ['wp_peak', 'wp_4060']:
            notetxt = 'Values filled where wp_peak, wp_4060 >= %f' % metadata['maximum_wp'] + '. '

            if 'note' in ds[var].attrs:
                ds[var].attrs['note'] = notetxt + ds[var].attrs['note']
            else:
                ds[var].attrs.update({'note': notetxt})

    return ds


def trim_min_wh(ds, metadata):
    """
    QA/QC
    Trim wave data based on minimum wave height as specified in metadata
    """

    if 'minimum_wh' in metadata:
        print('Trimming using minimum wave height of %f m'
            % metadata['minimum_wh'])
        ds = ds.where(ds['wh_4061'] > metadata['minimum_wh'])

        for var in ['wp_peak', 'wp_4060', 'wh_4061']:
            notetxt = 'Values filled where wh_4061 <= %f' % metadata['minimum_wh'] + '. '

            if 'note' in ds[var].attrs:
                ds[var].attrs['note'] = notetxt + ds[var].attrs['note']
            else:
                ds[var].attrs.update({'note': notetxt})

    return ds


def trim_wp_ratio(ds, metadata):
    """
    QA/QC
    Trim wave data based on maximum ratio of wp_peak to wp_4060
    """

    if 'wp_ratio' in metadata:
        print('Trimming using maximum ratio of wp_peak to wp_4060 of %f'
            % metadata['wp_ratio'])
        for var in ['wp_peak', 'wp_4060']:
            ds[var] = ds[var].where(ds['wp_peak']/ds['wp_4060'] < metadata['wp_ratio'])

        for var in ['wp_peak', 'wp_4060']:
            notetxt = 'Values filled where wp_peak:wp_4060 >= %f' % metadata['wp_ratio'] + '. '

            if 'note' in ds[var].attrs:
                ds[var].attrs['note'] = notetxt + ds[var].attrs['note']
            else:
                ds[var].attrs.update({'note': notetxt})

    return ds


def rename_time(nc_filename):
    """
    Rename time variables. Need to use netCDF4 module since xarray seems to have
    issues with the naming of time variables/dimensions
    """

    nc = netCDF4.Dataset(nc_filename, 'r+')
    timebak = nc['epic_time'][:]
    nc.renameVariable('time', 'time_cf')
    nc.renameVariable('epic_time', 'time')
    nc.renameVariable('epic_time2', 'time2')
    nc.close()

    # need to do this in two steps after renaming the variable
    # not sure why, but it works this way
    nc = netCDF4.Dataset(nc_filename, 'r+')
    nc['time'][:] = timebak
    nc.close()


def ds_add_attrs(ds, metadata):
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

    def add_attributes(var, metadata, INFO):
        var.attrs.update({'serial_number': INFO['serial_number'],
            'initial_instrument_height': metadata['initial_instrument_height'],
            # 'nominal_instrument_depth': metadata['nominal_instrument_depth'], # FIXME
            'height_depth_units': 'm',
            'sensor_type': INFO['INST_TYPE'],
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
        add_attributes(ds[var], metadata, ds.attrs)
        ds[var].attrs.update({'minimum': ds[var].min().values,
            'maximum': ds[var].max().values})

    return ds

def main():
    import argparse
    import yaml

    parser = argparse.ArgumentParser(description='Convert processed .nc files using DIWASP')
    parser.add_argument('gatts', help='path to global attributes file (gatts formatted)')
    parser.add_argument('config', help='path to ancillary config file (YAML formatted)')

    args = parser.parse_args()

    # initialize metadata from the globalatts file
    metadata = rsklib.read_globalatts(args.gatts)

    # Add additional metadata from metadata config file
    config = yaml.safe_load(open(args.config))

    for k in config:
        metadata[k] = config[k]

    ds = nc_to_diwasp(metadata)

    return ds

if __name__ == '__main__':
    main()
