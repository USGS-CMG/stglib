from __future__ import division, print_function

import xarray as xr

from .. import exo
from ..core import utils
from . import aqdutils


def cdf_to_nc(cdf_filename, atmpres=False):
    """
    Load a "raw" .cdf file and generate a processed .nc file
    """

    # Load raw .cdf data
    VEL = aqdutils.load_cdf(cdf_filename, atmpres=atmpres)

    # Clip data to in/out water times or via good_ens
    VEL = utils.clip_ds(VEL)

    # Create water_depth attribute
    # VEL = utils.create_water_depth(VEL)
    VEL = utils.create_nominal_instrument_depth(VEL)

    # Create depth variable depending on orientation
    VEL, T, T_orig = aqdutils.set_orientation(VEL, VEL["TransMatrix"].values)

    # Transform coordinates from, most likely, BEAM to ENU
    u, v, w = aqdutils.coord_transform(
        VEL["VEL1"].values,
        VEL["VEL2"].values,
        VEL["VEL3"].values,
        VEL["Heading"].values,
        VEL["Pitch"].values,
        VEL["Roll"].values,
        T,
        T_orig,
        VEL.attrs["AQDCoordinateSystem"],
    )

    VEL["U"] = xr.DataArray(u, dims=("time", "bindist"))
    VEL["V"] = xr.DataArray(v, dims=("time", "bindist"))
    VEL["W"] = xr.DataArray(w, dims=("time", "bindist"))

    VEL = aqdutils.magvar_correct(VEL)

    VEL["AGC"] = (VEL["AMP1"] + VEL["AMP2"] + VEL["AMP3"]) / 3

    VEL = aqdutils.trim_vel(VEL)

    VEL = aqdutils.make_bin_depth(VEL)

    # Reshape and associate dimensions with lat/lon
    for var in [
        "U",
        "V",
        "W",
        "AGC",
        "Pressure",
        "Temperature",
        "Heading",
        "Pitch",
        "Roll",
        "bin_depth",
        "Pressure_ac",
    ]:
        if var in VEL:
            VEL = utils.add_lat_lon(VEL, var)

    # swap_dims from bindist to depth
    VEL = ds_swap_dims(VEL)

    # Rename DataArrays for EPIC compliance
    VEL = aqdutils.ds_rename(VEL)

    # Drop unused variables
    VEL = ds_drop(VEL)

    # Add EPIC and CMG attributes
    VEL = aqdutils.ds_add_attrs(VEL)

    # should function this
    for var in VEL.data_vars:
        VEL = exo.trim_max_diff(VEL, var)
        VEL = exo.trim_min_diff(VEL, var)
        VEL = exo.trim_min(VEL, var)
        VEL = exo.trim_max(VEL, var)

    # Add min/max values
    VEL = utils.add_min_max(VEL)

    # Add DELTA_T for EPIC compliance
    VEL = aqdutils.add_delta_t(VEL)

    # Add start_time and stop_time attrs
    VEL = utils.add_start_stop_time(VEL)

    # Add history showing file used
    VEL = utils.add_history(VEL)

    # Rename time variables for EPIC compliance, keeping a time_cf coorindate.
    if utils.is_cf(VEL):
        pass
    else:
        VEL = utils.rename_time(VEL)

    if utils.is_cf(VEL):
        VEL = utils.add_standard_names(VEL)

    for var in VEL.variables:
        if (var not in VEL.coords) and ("time" not in var):
            # cast as float32
            VEL = utils.set_var_dtype(VEL, var)

    if "prefix" in VEL.attrs:
        nc_filename = VEL.attrs["prefix"] + VEL.attrs["filename"] + "-a.nc"
    else:
        nc_filename = VEL.attrs["filename"] + "-a.nc"

    if utils.is_cf(VEL):
        VEL.to_netcdf(nc_filename, encoding={"time": {"dtype": "i4"}})
        utils.check_compliance(nc_filename, checker_names=["cf:1.6"])
    else:
        VEL.to_netcdf(nc_filename, unlimited_dims=["time"])

    print("Done writing netCDF file", nc_filename)

    return VEL


# TODO: add analog input variables (OBS, NTU, etc)


def ds_swap_dims(ds):
    ds = ds.swap_dims({"bindist": "depth"})
    # need to swap dims and then reassign bindist to be a normal variable
    # (no longer a coordinate)
    valbak = ds["bindist"].values
    ds = ds.drop("bindist")
    ds["bindist"] = xr.DataArray(valbak, dims="depth")

    return ds


def ds_drop(ds):
    """
    Drop old DataArrays from Dataset that won't make it into the final .nc file
    """

    todrop = [
        "VEL1",
        "VEL2",
        "VEL3",
        "AMP1",
        "AMP2",
        "AMP3",
        "TransMatrix",
        "AnalogInput1",
        "AnalogInput2",
        "jd",
        "Depth",
    ]

    if ("AnalogInput1" in ds.attrs) and (ds.attrs["AnalogInput1"] == "true"):
        todrop.remove("AnalogInput1")

    if ("AnalogInput2" in ds.attrs) and (ds.attrs["AnalogInput2"] == "true"):
        todrop.remove("AnalogInput2")

    return ds.drop([t for t in todrop if t in ds.variables])
