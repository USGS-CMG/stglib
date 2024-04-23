import warnings

import numpy as np
import xarray as xr

from ..aqd import aqdutils
from ..core import qaqc, utils


def cdf_to_nc(cdf_filename, atmpres=False):
    """
    Load a raw .cdf file and generate a processed .nc file
    """

    # Load raw .cdf data
    ds = aqdutils.load_cdf(cdf_filename, atmpres=atmpres)

    # Clip data to in/out water times or via good_ens
    ds = utils.clip_ds(ds)

    ds = utils.create_nominal_instrument_depth(ds)

    ds, T, T_orig = set_orientation(ds, ds["TransMatrix"].values)

    u, v, w = aqdutils.coord_transform(
        ds["VEL1"],
        ds["VEL2"],
        ds["VEL3"],
        ds["Heading"].values,
        ds["Pitch"].values,
        ds["Roll"].values,
        T,
        T_orig,
        ds.attrs["VECCoordinateSystem"],
    )

    ds["U"] = xr.DataArray(u, dims=("time", "sample"))
    ds["V"] = xr.DataArray(v, dims=("time", "sample"))
    ds["W"] = xr.DataArray(w, dims=("time", "sample"))

    ds = aqdutils.magvar_correct(ds)

    # Rename DataArrays for EPIC compliance
    ds = aqdutils.ds_rename(ds)

    ds = dist_to_boundary(ds)

    ds = scale_analoginput(ds)

    # Drop unused variables
    ds = ds_drop(ds)

    ds = qaqc.drop_vars(ds)

    # Add EPIC and CMG attributes
    ds = aqdutils.ds_add_attrs(ds, inst_type="VEC")

    ds = associate_z_coord(ds)

    for v in ds.data_vars:
        # need to do this or else a "coordinates" attribute with value of "burst" hangs around
        ds[v].encoding["coordinates"] = None
        ds = qaqc.trim_warmup(ds, v)

    # Add start_time and stop_time attrs
    ds = utils.add_start_stop_time(ds)

    # Add history showing file used
    ds = utils.add_history(ds)

    ds = utils.add_standard_names(ds)

    if "prefix" in ds.attrs:
        nc_filename = ds.attrs["prefix"] + ds.attrs["filename"] + "-a.nc"
    else:
        nc_filename = ds.attrs["filename"] + "-a.nc"

    ds.to_netcdf(nc_filename, encoding={"time": {"dtype": "i4"}})
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])

    print("Done writing netCDF file", nc_filename)

    print("Creating burst-mean statistics file")

    dsmean = ds.mean(dim="sample", keep_attrs=True)

    # we already applied ClockDrift, ClockError in dat2cdf so don't re-apply it here.
    # but we do want to shift time since we are presenting mean values.
    dsmean = utils.shift_time(
        dsmean,
        dsmean.attrs["VECSamplesPerBurst"] / dsmean.attrs["VECSamplingRate"] / 2,
        apply_clock_error=False,
        apply_clock_drift=False,
    )

    if "prefix" in dsmean.attrs:
        nc_filename = dsmean.attrs["prefix"] + dsmean.attrs["filename"] + "-s.nc"
    else:
        nc_filename = dsmean.attrs["filename"] + "-s.nc"

    dsmean.to_netcdf(nc_filename, encoding={"time": {"dtype": "i4"}})
    utils.check_compliance(nc_filename, conventions=dsmean.attrs["Conventions"])

    print("Done writing netCDF file", nc_filename)

    return ds


def set_orientation(VEL, T):
    """
    Create z variable depending on instrument orientation
    Mostly taken from ../aqd/aqdutils.py
    """

    if "Pressure_ac" in VEL:
        presvar = "Pressure_ac"
    else:
        presvar = "Pressure"

    geopotential_datum_name = None

    if "NAVD88_ref" in VEL.attrs or "NAVD88_elevation_ref" in VEL.attrs:
        # if we have NAVD88 elevations of the bed, reference relative to the instrument height in NAVD88
        if "NAVD88_ref" in VEL.attrs:
            navd88_ref = VEL.attrs["NAVD88_ref"]
        elif "NAVD88_elevation_ref" in VEL.attrs:
            navd88_ref = VEL.attrs["NAVD88_elevation_ref"]

        # elev = VEL.attrs["NAVD88_ref"] + VEL.attrs["transducer_offset_from_bottom"]
        elev_vel = navd88_ref + VEL.attrs["velocity_sample_volume_height"]
        elev_pres = navd88_ref + VEL.attrs["pressure_sensor_height"]
        if "AnalogInput1_height" in VEL.attrs:
            elev_ai1 = navd88_ref + VEL.attrs["AnalogInput1_height"]
        if "AnalogInput2_height" in VEL.attrs:
            elev_ai2 = navd88_ref + VEL.attrs["AnalogInput2_height"]

        long_name = "height relative to NAVD88"
        geopotential_datum_name = "NAVD88"
    elif "height_above_geopotential_datum" in VEL.attrs:
        # elev = (
        #     VEL.attrs["height_above_geopotential_datum"]
        #     + VEL.attrs["transducer_offset_from_bottom"]
        # )
        elev_vel = (
            VEL.attrs["height_above_geopotential_datum"]
            + VEL.attrs["velocity_sample_volume_height"]
        )
        elev_pres = (
            VEL.attrs["height_above_geopotential_datum"]
            + VEL.attrs["pressure_sensor_height"]
        )
        if "AnalogInput1_height" in VEL.attrs:
            elev_ai1 = (
                VEL.attrs["height_above_geopotential_datum"]
                + VEL.attrs["AnalogInput1_height"]
            )
        if "AnalogInput2_height" in VEL.attrs:
            elev_ai2 = (
                VEL.attrs["height_above_geopotential_datum"]
                + VEL.attrs["AnalogInput2_height"]
            )

        long_name = f"height relative to {VEL.attrs['geopotential_datum_name']}"
        geopotential_datum_name = VEL.attrs["geopotential_datum_name"]
    else:
        # if we don't have NAVD88 elevations, reference to sea-bed elevation
        # elev = VEL.attrs["transducer_offset_from_bottom"]
        elev_vel = VEL.attrs["velocity_sample_volume_height"]
        elev_pres = VEL.attrs["pressure_sensor_height"]
        if "AnalogInput1_height" in VEL.attrs:
            elev_ai1 = VEL.attrs["AnalogInput1_height"]
        if "AnalogInput2_height" in VEL.attrs:
            elev_ai2 = VEL.attrs["AnalogInput2_height"]

        long_name = "height relative to sea bed"

    T_orig = T.copy()

    # User orientation refers to probe orientation
    # Nortek status code orientation refers to z-axis positive direction
    # See Nortek "The Comprehensive Manual - Velocimeters"
    # section 3.1.7 Orientation of Vector probes
    userorient = VEL.attrs["orientation"]
    # last bit of statuscode is orientation
    sc = str(VEL["StatusCode"].isel(time=int(len(VEL["time"]) / 2)).values)[-1]
    if sc == "0":
        scname = "UP"
    elif sc == "1":
        scname = "DOWN"
    headtype = VEL.attrs["VECHeadSerialNumber"][0:3]

    print(
        f"Instrument reported {headtype} case with orientation status code {sc} -> z-axis positive {scname} at middle of deployment"
    )

    if userorient == "UP":
        print("User instructed probe is pointing UP (sample volume above probe)")
    elif userorient == "DOWN":
        print("User instructed probe is pointing DOWN (sample volume below probe)")
    else:
        raise ValueError("Could not determine instrument orientation from user input")

    flag = False
    if headtype == "VEC":
        if sc == "0" and userorient == "UP":
            flag = True
        elif sc == "1" and userorient == "DOWN":
            flag = True
    elif headtype == "VCH":
        if sc == "0" and userorient == "DOWN":
            flag = True
        elif sc == "1" and userorient == "UP":
            flag = True

    if flag is False:
        print(
            "User-provided orientation matches orientation status code at middle of deployment"
        )
    elif flag is True:
        warnings.warn(
            "User-provided orientation does not match orientation status code at middle of deployment"
        )

        histtext = "Modifying transformation matrix to match user-provided orientation"

        warnings.warn(histtext)

        VEL = utils.insert_history(VEL, histtext)

        T[1, :] = -T[1, :]
        T[2, :] = -T[2, :]

    diff = elev_pres - elev_vel
    VEL["depthvel"] = xr.DataArray(np.nanmean(VEL[presvar]) + [diff], dims="depthvel")
    VEL["depthpres"] = xr.DataArray([np.nanmean(VEL[presvar])], dims="depthpres")
    VEL["zvel"] = xr.DataArray([elev_vel], dims="zvel")
    VEL["zpres"] = xr.DataArray([elev_pres], dims="zpres")
    if "AnalogInput1_height" in VEL.attrs:
        VEL["zai1"] = xr.DataArray([elev_ai1], dims="zai1")
    if "AnalogInput2_height" in VEL.attrs:
        VEL["zai2"] = xr.DataArray([elev_ai2], dims="zai2")

    lnshim = {
        "zvel": "of velocity sensor",
        "zpres": "of pressure sensor",
        "zai1": "of analog input 1",
        "zai2": "of analog input 2",
    }
    for z in ["zvel", "zpres", "zai1", "zai2"]:
        if z not in VEL:
            continue
        VEL[z].attrs["standard_name"] = "height"
        VEL[z].attrs["units"] = "m"
        VEL[z].attrs["positive"] = "up"
        VEL[z].attrs["axis"] = "Z"
        VEL[z].attrs["long_name"] = f"{long_name} {lnshim[z]}"
        if geopotential_datum_name:
            VEL[z].attrs["geopotential_datum_name"] = geopotential_datum_name

    for d in ["depthvel", "depthpres"]:
        VEL[d].attrs["standard_name"] = "depth"
        VEL[d].attrs["units"] = "m"
        VEL[d].attrs["positive"] = "down"
        VEL[d].attrs[
            "long_name"
        ] = f"depth {lnshim[z]} below mean sea level of deployment"

    # "z" is ambiguous, so drop it from Dataset for now
    # FIXME: remove creation of z variable above instead of just dropping it
    # VEL = VEL.drop("z")

    return VEL, T, T_orig


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
        "Depth",
        "Checksum",
        "ErrorCode",
        "StatusCode",
        "vel1_1277",
        "vel2_1278",
        "vel3_1279",
    ]

    if ("AnalogInput1" in ds.attrs) and (ds.attrs["AnalogInput1"].lower() == "true"):
        todrop.remove("AnalogInput1")

    if ("AnalogInput2" in ds.attrs) and (ds.attrs["AnalogInput2"].lower() == "true"):
        todrop.remove("AnalogInput2")

    return ds.drop([t for t in todrop if t in ds.variables])


def scale_analoginput(ds):
    """convert AnalogInput from counts to volts"""
    ds["AnalogInput1"] = ds["AnalogInput1"] * 5 / 65535
    notetxt = "Converted from counts to volts: volts=counts*5/65535."
    ds = utils.insert_note(ds, "AnalogInput1", notetxt)
    ds["AnalogInput2"] = ds["AnalogInput2"] * 5 / 65535
    notetxt = "Converted from counts to volts: volts=counts*5/65535."
    ds = utils.insert_note(ds, "AnalogInput2", notetxt)

    return ds


def associate_z_coord(ds):
    """Associate the appropriate z coordinate to data variables.
    We do this because there are multiple relevant elevations per deployment
    (e.g., velocity and pressure were collected at different elevations,
    and we need to indicate this)"""

    for v in [
        "u_1205",
        "v_1206",
        "w_1204",
        "AGC1_1221",
        "AGC2_1222",
        "AGC3_1223",
        "SNR1",
        "SNR2",
        "SNR3",
        "cor1_1285",
        "cor2_1286",
        "cor3_1287",
    ]:
        if v in ds:
            # pass axis=-1 to add z dim to end for CF compliance
            ds[v] = ds[v].expand_dims("zvel", axis=-1)

    for v in ["P_1ac", "P_1"]:
        if v in ds:
            ds[v] = ds[v].expand_dims("zpres", axis=-1)

    if "AnalogInput1" in ds:
        ds["AnalogInput1"] = ds["AnalogInput1"].expand_dims("zai1", axis=-1)

    if "AnalogInput2" in ds:
        ds["AnalogInput2"] = ds["AnalogInput2"].expand_dims("zai2", axis=-1)

    return ds


def dist_to_boundary(ds):
    """Create range to boundary variable from start/end values"""
    ds["brange"] = (ds["DistProbeStartAvg"] + ds["DistProbeEndAvg"]) / 2
    ds["vrange"] = (ds["DistSVolStartAvg"] + ds["DistSVolEndAvg"]) / 2

    for v in [
        "DistProbeStartAvg",
        "DistProbeEndAvg",
        "DistSVolStartAvg",
        "DistSVolEndAvg",
    ]:
        ds = ds.drop(v)

    return ds
