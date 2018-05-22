from __future__ import division, print_function
import xarray as xr
from ..core import utils, waves


def nc_to_waves(nc_filename):

    ds = xr.open_dataset(nc_filename, autoclose=True, decode_times=False)

    ds = utils.epic_to_cf_time(ds)

    ds = utils.create_epic_time(ds)

    spec = waves.make_waves_ds(ds)

    for k in ['wp_peak', 'wh_4061', 'wp_4060', 'pspec']:
        ds[k] = spec[k]

    ds = utils.create_water_depth(ds)

    ds = ds.drop(['P_1', 'P_1ac', 'sample', 'Depth'])

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

    ds.to_netcdf(nc_filename, unlimited_dims=['time'])

    return ds
