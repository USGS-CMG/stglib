#!/usr/bin/env python

from __future__ import division, print_function
import xarray as xr
import sys
sys.path.append('/Users/dnowacki/Documents/rsklib')
import rsklib


def nc_to_diwasp(metadata):

    ds = xr.open_dataset(metadata['filename'] + '-wvsb-cal.nc', autoclose=True, decode_times=False)
    print(ds)
    ds['time'] = ds['time_cf']
    ds = ds.drop(['time_cf', 'time2'])
    ds = xr.decode_cf(ds, decode_times=True)

    ds = rsklib.rskcdf2nc.create_epic_time(ds)

    mat = xr.open_dataset(metadata['filename'] + '-wvs-diwasp-pres.nc', autoclose=True)

    for k in ['wp_peak', 'wh_4061', 'wp_4060']:
        ds[k] = xr.DataArray(mat[k], dims='time')

    ds['frequency'] = xr.DataArray(mat['frequency'], dims=('frequency'))

    ds['pspec'] = xr.DataArray(mat['pspec'], dims=('time', 'frequency'))

    ds, metadata = rsklib.rsknc2diwasp.create_water_depth(ds, metadata)

    ds = ds.drop(['P_1', 'P_1ac', 'sample'])

    ds = rsklib.rsknc2diwasp.trim_max_wp(ds, metadata)

    ds = rsklib.rsknc2diwasp.trim_min_wh(ds, metadata)

    ds = rsklib.rsknc2diwasp.trim_wp_ratio(ds, metadata)

    # Add attrs
    ds = rsklib.rsknc2diwasp.ds_add_attrs(ds, metadata)

    ds = rsklib.write_metadata(ds, metadata)

    write_nc(ds, metadata)

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
