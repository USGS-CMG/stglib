from __future__ import division, print_function
import sys
import netCDF4
import xarray as xr
from ..core import utils

def nc_to_diwasp(nc_filename):

    ds = xr.open_dataset(nc_filename, autoclose=True, decode_times=False)

    ds = utils.epic_to_cf_time(ds)

    ds = utils.create_epic_time(ds)

    mat = xr.open_dataset(ds.attrs['filename'][:-2] + 'diwasp.nc', autoclose=True)

    for k in ['wp_peak', 'wh_4061', 'wp_4060']:
        ds[k] = xr.DataArray(mat[k], dims='time')

    ds['frequency'] = xr.DataArray(mat['frequency'], dims=('frequency'))

    ds['pspec'] = xr.DataArray(mat['pspec'], dims=('time', 'frequency'))

    ds = utils.create_water_depth(ds)

    ds = ds.drop(['P_1', 'P_1ac', 'sample'])

    ds = utils.trim_max_wp(ds)

    ds = utils.trim_min_wh(ds)

    ds = utils.trim_max_wh(ds)

    ds = utils.trim_wp_ratio(ds)

    # Add attrs
    ds = utils.ds_add_attrs(ds)

    # Reshape and associate dimensions with lat/lon
    for var in ['wp_peak', 'wh_4061', 'wp_4060', 'pspec']:
        if var in ds:
            ds = utils.add_lat_lon(ds, var)

    # assign min/max (need to do this after trimming):
    ds = utils.add_min_max(ds)

    ds = utils.ds_add_diwasp_history(ds)

    nc_filename = ds.attrs['filename'] + 's-a.nc'

    ds = utils.rename_time(ds)

    ds.to_netcdf(nc_filename)

    return ds
