from __future__ import division, print_function

import numpy as np
import xarray as xr

from ..core import utils, waves


def nc_to_diwasp(nc_filename, format='NETCDF3_64BIT'):

    ds = utils.open_time_2d_dataset(nc_filename)

    ds = utils.epic_to_cf_time(ds)

    ds = utils.create_epic_times(ds)

    mat = xr.open_dataset(ds.attrs['filename'] + 'wvs-diwasp.nc').load()
    mat.close()

    ds['frequency'] = xr.DataArray(mat['frequency'], dims=('frequency'))

    # Note we convert the polar/cartesian coordinates provided by DIWASP to
    # compass coordinates
    dirs = waves.polar2compass(waves.to2from(mat['direction']))

    # Return only unique directions (diwasp has double directions...)
    # FIXME: Why is this?
    _, idx = np.unique(dirs, return_index=True)

    ds['direction'] = xr.DataArray(dirs[idx], dims=('direction'))
    ds['direction'].encoding['_FillValue'] = None

    ds['dspec'] = xr.DataArray(mat['dspec'][:, idx, :],
                               dims=('time', 'direction', 'frequency'))

    pspec = np.trapz(ds['dspec'].values, x=ds['direction'], axis=1)

    m0 = waves.make_moment(ds['frequency'].values, pspec, 0)
    m2 = waves.make_moment(ds['frequency'].values, pspec, 2)

    ds['wh_4061'] = xr.DataArray(waves.make_Hs(m0), dims='time')

    ds['wp_4060'] = xr.DataArray(waves.make_Tm(m0, m2), dims='time')

    for k in ['wp_peak']:
        ds[k] = xr.DataArray(mat[k], dims='time')

    for k in ['wvdir', 'dwvdir']:
        ds[k] = xr.DataArray(waves.polar2compass(waves.to2from(mat[k])),
                             dims='time')

    ds['pspec'] = xr.DataArray(pspec, dims=('time', 'frequency'))

    ds['wd_4062'] = xr.DataArray(waves.polar2compass(waves.to2from(
                                     waves.make_mwd(ds['frequency'].values,
                                                    ds['direction'].values,
                                                    ds['dspec'].values))),
                                 dims='time')

    # ds = utils.create_water_depth(ds)

    ds = utils.create_water_depth_var(ds)

    # Remove old variables as we just want to keep the wave statistics
    for v in ['P_1',
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
              'soundspeed',
              'Battery',
              'Hdg_1215',
              'Ptch_1216',
              'Roll_1217',
              'Bat_106',
              'Depth']:
        if v in ds:
            ds = ds.drop(v)

    # remove depth from bin_depth
    ds['bin_depth'] = ds['bin_depth'].squeeze(dim='depth')

    ds = utils.trim_max_wp(ds)

    ds = utils.trim_min_wh(ds)

    ds = utils.trim_max_wh(ds)

    ds = utils.trim_wp_ratio(ds)

    # Add attrs
    ds = utils.ds_add_attrs(ds)

    ds = utils.ds_add_diwasp_history(ds)

    # add lat/lon coordinates to each variable
    for var in ds.variables:
        if (var not in ds.coords) and ('time' not in var):
            ds = utils.add_lat_lon(ds, var)
            # cast as float32
            ds = utils.set_var_dtype(ds, var)

    # remove lat and lon from burst dimension: they are singleton so can
    # remove them both with one call to squeeze()
    ds['burst'] = ds['burst'].squeeze()

    nc_filename = ds.attrs['filename'] + 'wvs_diwasp-cal.nc'

    ds = utils.rename_time(ds)

    print('Writing to netCDF')

    ds.to_netcdf(nc_filename, format=format)

    print('Done creating', nc_filename)

    return ds
