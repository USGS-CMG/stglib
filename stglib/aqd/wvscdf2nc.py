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

    ds = qaqc.make_bin_depth(ds)

    ds = qaqc.ds_rename(ds, waves=True)

    ds = qaqc.ds_add_attrs(ds, waves=True)

    ds = utils.add_min_max(ds)

    nc_filename = ds.attrs['filename'] + 'wvsb-cal.nc'

    ds.to_netcdf(nc_filename, unlimited_dims='time')
    print('Done writing netCDF file', nc_filename)

    # rename time variables after the fact to conform with EPIC/CMG standards
    utils.rename_time(nc_filename)

    return ds
