#!/usr/bin/env python

from __future__ import division, print_function
import xarray as xr
import sys
from ..core import utils


def nc_to_diwasp(nc_filename):

    ds = xr.open_dataset(nc_filename, autoclose=True, decode_times=False)

    ds = utils.epic_to_cf_time(ds)

    ds = utils.create_epic_time(ds)

    mat = xr.open_dataset(ds.attrs['filename'] + 'wvs-diwasp.nc', autoclose=True)

    ds['frequency'] = xr.DataArray(mat['frequency'], dims=('frequency'))

    ds['direction'] = xr.DataArray(mat['direction'], dims=('direction'))

    for k in ['wp_peak', 'wh_4061', 'wp_4060', 'wvdir', 'dwvdir']:
        ds[k] = xr.DataArray(mat[k], dims='time')

    ds['pspec'] = xr.DataArray(mat['pspec'], dims=('time', 'frequency'))

    ds['dspec'] = xr.DataArray(mat['dspec'], dims=('time', 'direction', 'frequency'))

    ds = create_water_depth(ds)

    # Remove old variables as we just want to keep the wave statistics
    ds = ds.drop(['P_1',
                  'P_1ac',
                  'sample',
                  'Tx_1211',
                  'vel1_1277',
                  'vel2_1278',
                  'vel3_1279',
                  'U',
                  'V',
                  'W',
                  'avgamp1',
                  'avgamp2',
                  'avgamp3',
                  'AGC1_1221',
                  'AGC2_1222',
                  'AGC3_1223',
                  'TransMatrix',
                  'nrecs',
                  'burst',
                  'soundspeed',
                  'Battery',
                  'Hdg_1215',
                  'Ptch_1216',
                  'Roll_1217'])

    ds = utils.trim_max_wp(ds)

    ds = utils.trim_min_wh(ds)

    ds = utils.trim_max_wh(ds)

    ds = utils.trim_wp_ratio(ds)

    # Add attrs
    ds = utils.ds_add_attrs(ds)

    ds = utils.ds_add_diwasp_history(ds)

    nc_filename = ds.attrs['filename'] + 'wvs-a.nc'

    ds = utils.rename_time(ds)

    ds.to_netcdf(nc_filename)

    print('Done creating', nc_filename)

    return ds


def create_water_depth(ds):
    """Create water_depth variable"""

    if 'initial_instrument_height' in ds.attrs:
        if 'P_1ac' in ds:
            ds.attrs['nominal_instrument_depth'] = ds['P_1ac'].mean(dim='sample').values
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


def write_nc(ds, metadata):
    """Write cleaned and trimmed Dataset to .nc file"""

    nc_filename = metadata['filename'] + 'wvs_diwasp-cal.nc'


    ds.to_netcdf(nc_filename, unlimited_dims='time', engine='netcdf4')

    # rename time variables after the fact to conform with EPIC/CMG standards
    rename_time(nc_filename)

def main():
    import sys
    sys.path.insert(0, '/Users/dnowacki/Documents/aqdlib')
    import aqdlib
    import argparse
    import yaml


    parser = argparse.ArgumentParser(description='Convert processed .nc files using DIWASP')
    parser.add_argument('gatts', help='path to global attributes file (gatts formatted)')
    parser.add_argument('config', help='path to ancillary config file (YAML formatted)')

    args = parser.parse_args()

    # initialize metadata from the globalatts file
    metadata = aqdlib.read_globalatts(args.gatts)

    # Add additional metadata from metadata config file
    config = yaml.safe_load(open(args.config))

    for k in config:
        metadata[k] = config[k]

    ds = nc_to_diwasp(metadata)

if __name__ == '__main__':
    main()
