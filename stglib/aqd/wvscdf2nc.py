from __future__ import division, print_function

import xarray as xr
from ..core import utils
from . import qaqc


def cdf_to_nc(cdf_filename, atmpres=False, writefile=True):

    # Load raw .cdf data
    ds = qaqc.load_cdf(cdf_filename, atmpres=atmpres)

    # Clip data to in/out water times or via good_ens
    ds = utils.clip_ds(ds)

    # Create water_depth variables
    ds = utils.create_water_depth(ds)

    # Create depth variable depending on orientation
    ds, T, T_orig = qaqc.set_orientation(ds, ds['TransMatrix'].values)

    # Transform coordinates from, most likely, BEAM to ENU
    u, v, w = qaqc.coord_transform(ds['VEL1'].values,
                                   ds['VEL2'].values,
                                   ds['VEL3'].values,
                                   ds['Heading'].values,
                                   ds['Pitch'].values,
                                   ds['Roll'].values,
                                   T,
                                   T_orig,
                                   ds.attrs['AQDCoordinateSystem'])

    ds['U'] = xr.DataArray(u, dims=('time', 'sample'))
    ds['V'] = xr.DataArray(v, dims=('time', 'sample'))
    ds['W'] = xr.DataArray(w, dims=('time', 'sample'))

    ds = qaqc.magvar_correct(ds)

    ds = qaqc.make_bin_depth(ds)

    ds = qaqc.ds_rename(ds, waves=True)

    ds = qaqc.ds_add_attrs(ds, waves=True)

    ds = utils.add_min_max(ds)

    for var in ds.variables:
        if (var not in ds.coords) and ('time' not in var):
            # cast as float32
            ds = utils.set_var_dtype(ds, var)

    ds = utils.ds_coord_no_fillvalue(ds)

    if writefile:
        nc_filename = ds.attrs['filename'] + 'wvsb-cal.nc'
        ds.to_netcdf(nc_filename)

        # Rename time variables for EPIC compliance, keeping a time_cf
        # coorindate.
        utils.rename_time_2d(nc_filename)

        print('Done writing netCDF file', nc_filename)

    return ds
