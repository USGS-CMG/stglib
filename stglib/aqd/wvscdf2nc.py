from __future__ import division, print_function

import xarray as xr
from ..core import utils
from . import qaqc

def cdf_to_nc(cdf_filename, atmpres=False):

    # Load raw .cdf data
    VEL = qaqc.load_cdf(cdf_filename, atmpres=atmpres)

    # Clip data to in/out water times or via good_ens
    VEL = utils.clip_ds(VEL)

    # Create water_depth variables
    VEL = utils.create_water_depth(VEL)

    # Create depth variable depending on orientation
    VEL, T = qaqc.set_orientation(VEL, VEL['TransMatrix'].values)

    VEL = qaqc.make_bin_depth(VEL)

    VEL = qaqc.ds_rename(VEL, waves=True)

    VEL = qaqc.ds_add_attrs(VEL, waves=True)

    VEL = utils.add_min_max(VEL)

    nc_filename = VEL.attrs['filename'] + 'wvsb-cal.nc'

    VEL.to_netcdf(nc_filename, unlimited_dims='time')
    print('Done writing netCDF file', nc_filename)

    # rename time variables after the fact to conform with EPIC/CMG standards
    utils.rename_time(nc_filename)

    return VEL
