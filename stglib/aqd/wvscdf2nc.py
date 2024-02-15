from ..core import utils
from . import aqdutils, cdf2nc


def cdf_to_nc(cdf_filename, atmpres=False, writefile=True):
    """
    Load a raw .cdf file and generate a processed .nc file
    """

    # Load raw .cdf data
    ds = aqdutils.load_cdf(cdf_filename, atmpres=atmpres)

    # Clip data to in/out water times or via good_ens
    ds = utils.clip_ds(ds, wvs=True)

    # Create water_depth variables
    # ds = utils.create_water_depth(ds)
    ds = utils.create_nominal_instrument_depth(ds)

    # Create depth variable depending on orientation
    ds, T, T_orig = aqdutils.set_orientation(ds, ds["TransMatrix"].values)

    ds = aqdutils.apply_wave_coord_output(ds, T, T_orig)

    # Make bin_depth variable
    ds = aqdutils.make_bin_depth(ds, waves=True)

    ds = cdf2nc.ds_swap_dims(ds)
    # Rename DataArrays within Dataset for EPIC compliance
    # and append depth coord to velocities and amplitudes
    ds = aqdutils.ds_rename(ds, waves=True)
    # add EPIC and CMG attributes, set _FillValue
    ds = aqdutils.ds_add_attrs(ds, waves=True)

    # Add minimum and maximum attributes
    ds = utils.add_min_max(ds)

    # Add DELTA_T for EPIC compliance
    ds = aqdutils.add_delta_t(ds, waves=True)

    # Add start_time and stop_time attrs
    ds = utils.add_start_stop_time(ds)

    # Add history showing file used
    ds = utils.add_history(ds)

    ds = utils.add_standard_names(ds)

    # Cast vars as float32
    for var in ds.variables:
        if (var not in ds.coords) and ("time" not in var):
            # cast as float32
            ds = utils.set_var_dtype(ds, var)

    # Ensure no _FillValue is assigned to coordinates
    ds = utils.ds_coord_no_fillvalue(ds)

    if writefile:
        nc_filename = ds.attrs["filename"] + "wvb-cal.nc"
        ds.to_netcdf(nc_filename, encoding={"time": {"dtype": "i4"}})
        utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])
        print("Done writing netCDF file", nc_filename)

    return ds
