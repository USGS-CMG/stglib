from __future__ import division, print_function
import xarray as xr
from ..core import utils, waves


def nc_to_waves(nc_filename):

    ds = xr.load_dataset(nc_filename, decode_times=False)

    if utils.is_cf(ds):
        for k in ds:
            if '_time' in k:
                ds = ds.drop(k)
        ds = xr.decode_cf(ds)
    else:
        ds = utils.epic_to_cf_time(ds)
        ds = utils.create_epic_times(ds)

    spec = waves.make_waves_ds(ds)

    for k in ['wp_peak', 'wh_4061', 'wp_4060', 'pspec']:
        ds[k] = spec[k]

    # ds = utils.create_water_depth(ds)

    # Remove old variables as we just want to keep the wave statistics
    keys = ['P_1',
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
             'Roll_1217']

    for k in keys:
        if k in ds:
            ds = ds.drop(k)

    ds = utils.trim_max_wp(ds)

    ds = utils.trim_min_wh(ds)

    ds = utils.trim_max_wh(ds)

    ds = utils.trim_wp_ratio(ds)

    # Add attrs
    ds = utils.ds_add_attrs(ds)

    nc_filename = ds.attrs['filename'] + 'wvs-a.nc'

    ds.to_netcdf(nc_filename, unlimited_dims=['time'])

    print('Done creating', nc_filename)

    return ds
