from __future__ import division, print_function

import xarray as xr
from ..core import utils
from . import qaqc

def cdf_to_nc(cdf_filename, atmpres=False):

    # Load raw .cdf data
    ds = qaqc.load_cdf(cdf_filename, atmpres=atmpres)

    # Clip data to in/out water times or via good_ens
    ds = utils.clip_ds(ds)

    # Create water_depth variables
    ds = utils.create_water_depth(ds)

    # Create depth variable depending on orientation
    ds, T = qaqc.set_orientation(ds, ds['TransMatrix'].values)

    # Transform coordinates from, most likely, BEAM to ENU
    u, v, w = qaqc.coord_transform(ds['VEL1'].values, ds['VEL2'].values, ds['VEL3'].values,
        ds['Heading'].values, ds['Pitch'].values, ds['Roll'].values, T, ds.attrs['AQDCoordinateSystem'])

    ds['U'] = xr.DataArray(u, dims=('time', 'sample'))
    ds['V'] = xr.DataArray(v, dims=('time', 'sample'))
    ds['W'] = xr.DataArray(w, dims=('time', 'sample'))

    ds = qaqc.magvar_correct(ds)

    ds = qaqc.make_bin_depth(ds)

    ds = qaqc.ds_rename(ds, waves=True)

    ds = qaqc.ds_add_attrs(ds, waves=True)

    ds = utils.add_min_max(ds)

    # Rename time variables for EPIC compliance, keeping a time_cf coorindate.
    ds = utils.rename_time(ds)

    nc_filename = ds.attrs['filename'] + 'wvsb-cal.nc'

    ds.to_netcdf(nc_filename, unlimited_dims='time')
    print('Done writing netCDF file', nc_filename)

    return ds
