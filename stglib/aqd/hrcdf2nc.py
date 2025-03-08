import xarray as xr

from ..core import filter, qaqc, utils
from . import aqdutils


def cdf_to_nc(cdf_filename, atmpres=False):
    """
    Load a raw .cdf file and generate a processed .nc file
    """

    # Load raw .cdf data
    VEL = aqdutils.load_cdf(cdf_filename, atmpres=atmpres)

    # Clip data to in/out water times or via good_ens
    VEL = utils.clip_ds(VEL)

    # Create water_depth attribute
    # VEL = utils.create_water_depth(VEL)
    VEL = utils.create_nominal_instrument_depth(VEL)

    # create Z depending on orientation
    VEL, T, T_orig = aqdutils.set_orientation(VEL, VEL["TransMatrix"].values)

    # Transform coordinates from, most likely, BEAM to ENU

    if "VEL1" in VEL:
        theu = VEL["VEL1"]
        thev = VEL["VEL2"]
        thew = VEL["VEL3"]
    elif "U" in VEL:
        theu = VEL["U"]
        thev = VEL["V"]
        thew = VEL["W"]
    elif "X" in VEL:
        theu = VEL["X"]
        thev = VEL["Y"]
        thew = VEL["Z"]

    out = "ENU"

    if VEL.attrs["AQDHRCoordinateSystem"] != out:
        histtext = f"Transforming data from {VEL.attrs['AQDHRCoordinateSystem']} coordinates to {out} coordinates."

        VEL = utils.insert_history(VEL, histtext)

    VEL["U"] = xr.zeros_like(theu)
    VEL["V"] = xr.zeros_like(thev)
    VEL["W"] = xr.zeros_like(thew)

    for n in range(len(VEL.sample)):
        u, v, w = aqdutils.coord_transform(
            theu.isel(sample=n).values,
            thev.isel(sample=n).values,
            thew.isel(sample=n).values,
            VEL["Heading"].isel(sample=n).values,
            VEL["Pitch"].isel(sample=n).values,
            VEL["Roll"].isel(sample=n).values,
            T,
            T_orig,
            VEL.attrs["AQDHRCoordinateSystem"],
            out=out,
        )
        VEL["U"][:, n, :] = u
        VEL["V"][:, n, :] = v
        VEL["W"][:, n, :] = w

    VEL = aqdutils.magvar_correct(VEL)

    VEL["AGC"] = (VEL["AMP1"] + VEL["AMP2"] + VEL["AMP3"]) / 3

    if "COR1" in VEL:
        VEL["COR"] = (VEL["COR1"] + VEL["COR2"] + VEL["COR3"]) / 3

    VEL = aqdutils.trim_vel(VEL)

    VEL = aqdutils.make_bin_depth(VEL)

    # swap vert dim to z or user specified in vert_dim
    VEL = aqdutils.ds_swap_dims(VEL)

    # Rename DataArrays for EPIC compliance
    VEL = aqdutils.ds_rename(VEL)

    # Drop unused variables
    VEL = ds_drop(VEL)

    VEL = qaqc.drop_vars(VEL)

    # Add EPIC and CMG attributes
    VEL = aqdutils.ds_add_attrs(VEL, hr=True)

    # should function this
    for var in VEL.data_vars:
        # check if any filtering before other qaqc
        VEL = filter.apply_butter_filt(VEL, var)
        VEL = filter.apply_med_filt(VEL, var)

        VEL = qaqc.trim_min(VEL, var)
        VEL = qaqc.trim_max(VEL, var)
        VEL = qaqc.trim_min_diff(VEL, var)
        VEL = qaqc.trim_min_diff_pct(VEL, var)
        VEL = qaqc.trim_max_diff(VEL, var)
        VEL = qaqc.trim_max_diff_pct(VEL, var)
        VEL = qaqc.trim_med_diff(VEL, var)
        VEL = qaqc.trim_med_diff_pct(VEL, var)
        VEL = qaqc.trim_bad_ens(VEL, var)
        VEL = qaqc.trim_maxabs_diff_2d(VEL, var)
        VEL = aqdutils.trim_single_bins(VEL, var)
        VEL = qaqc.trim_fliers(VEL, var)

    # after check for masking vars by other vars
    for var in VEL.data_vars:
        VEL = qaqc.trim_mask(VEL, var)
        VEL = qaqc.trim_mask_expr(VEL, var)

    # fill with AGC and Cor threshold
    VEL = aqdutils.fill_agc(VEL)
    VEL = aqdutils.fill_cor(VEL)

    # Add min/max values
    VEL = utils.add_min_max(VEL)

    if False:
        # Add DELTA_T for EPIC compliance
        VEL = aqdutils.add_delta_t(VEL)

    # Add start_time and stop_time attrs
    VEL = utils.add_start_stop_time(VEL)

    # Add history showing file used
    VEL = utils.add_history(VEL)

    VEL = utils.add_standard_names(VEL)

    if "prefix" in VEL.attrs:
        nc_filename = VEL.attrs["prefix"] + VEL.attrs["filename"] + "b-cal.nc"
    else:
        nc_filename = VEL.attrs["filename"] + "b-cal.nc"

    VEL.to_netcdf(nc_filename, encoding={"time": {"dtype": "i4"}})
    utils.check_compliance(nc_filename, conventions=VEL.attrs["Conventions"])

    print("Done writing netCDF file", nc_filename)

    # make burst averaged -a.nc file
    VEL = aqdutils.average_burst(VEL)

    # drop beam data
    for k in [
        "AGC1_1221",
        "AGC2_1222",
        "AGC3_1223",
        "cor1_1285",
        "cor2_1286",
        "cor3_1287",
        "vel1_1277",
        "vel2_1278",
        "vel3_1279",
    ]:
        if k in VEL:
            VEL = VEL.drop_vars(k)

    # assign min/max
    VEL = utils.add_min_max(VEL)

    VEL = utils.ds_coord_no_fillvalue(VEL)

    if "prefix" in VEL.attrs:
        nc_filename = VEL.attrs["prefix"] + VEL.attrs["filename"] + "-a.nc"
    else:
        nc_filename = VEL.attrs["filename"] + "-a.nc"

    VEL.to_netcdf(nc_filename, encoding={"time": {"dtype": "i4"}})
    utils.check_compliance(nc_filename, conventions=VEL.attrs["Conventions"])

    print("Done writing burst averaged netCDF file", nc_filename)

    return VEL


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
        "Burst",
        "Ensemble",
    ]

    if ("AnalogInput1" in ds.attrs) and (ds.attrs["AnalogInput1"].lower() == "true"):
        todrop.remove("AnalogInput1")

    if ("AnalogInput2" in ds.attrs) and (ds.attrs["AnalogInput2"].lower() == "true"):
        todrop.remove("AnalogInput2")

    return ds.drop([t for t in todrop if t in ds.variables])
