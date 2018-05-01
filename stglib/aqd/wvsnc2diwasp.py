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

    ds = utils.create_water_depth(ds)

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
    

def main():

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
