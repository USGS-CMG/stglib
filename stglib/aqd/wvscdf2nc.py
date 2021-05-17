from __future__ import division, print_function

from ..core import utils
from . import qaqc


def cdf_to_nc(cdf_filename, atmpres=False, writefile=True): # , format="NETCDF3_64BIT"): # don't think we need to fall back to netcdf3 any more

    # Load raw .cdf data
    ds = qaqc.load_cdf(cdf_filename, atmpres=atmpres)

    # Clip data to in/out water times or via good_ens
    ds = utils.clip_ds(ds, wvs=True)

    # Create water_depth variables
    # ds = utils.create_water_depth(ds)
    ds = utils.create_nominal_instrument_depth(ds)

    # Create depth variable depending on orientation
    ds, T, T_orig = qaqc.set_orientation(ds, ds["TransMatrix"].values)

    # Make bin_depth variable
    ds = qaqc.make_bin_depth(ds, waves=True)

    # Swap dimensions from bindist to depth
    ds = qaqc.swap_bindist_to_depth(ds)
    # Rename DataArrays within Dataset for EPIC compliance
    # and append depth coord to velocities and amplitudes
    ds = qaqc.ds_rename(ds, waves=True)

    # add EPIC and CMG attributes, set _FillValue
    ds = qaqc.ds_add_attrs(ds, waves=True)

    # Add DELTA_T for EPIC compliance
    ds = qaqc.add_delta_t(ds, waves=True)

    # Add minimum and maximum attributes
    ds = utils.add_min_max(ds)

    # Cast vars as float32
    for var in ds.variables:
        if (var not in ds.coords) and ("time" not in var):
            # cast as float32
            ds = utils.set_var_dtype(ds, var)

    # Need to add lat lon to certain variables
    for var in ["Hdg_1215", "Ptch_1216", "Roll_1217"]:
        ds = utils.add_lat_lon(ds, var)

    # Ensure no _FillValue is assigned to coordinates
    ds = utils.ds_coord_no_fillvalue(ds)

    if writefile:
        nc_filename = ds.attrs["filename"] + "wvsb-cal.nc"
        ds.to_netcdf(nc_filename)
        # Rename time variables for EPIC compliance, keeping a time_cf
        # coorindate.
        utils.rename_time_2d(nc_filename, ds)

        print("Done writing netCDF file", nc_filename)

    return ds
