import warnings

import numpy as np
import pandas as pd
import xarray as xr
from dask.diagnostics import ProgressBar

from ..aqd import aqdutils
from ..core import qaqc, transform, utils


def cdf_to_nc(cdf_filename, atmpres=False):
    """
    Load a raw .cdf file and generate a processed .nc file
    """

    # Load raw .cdf data
    ds = aqdutils.load_cdf(cdf_filename, atmpres=atmpres)

    # Clip data to in/out water times or via good_ens
    ds = utils.clip_ds(ds)

    ds = utils.create_nominal_instrument_depth(ds)

    ds = set_orientation(ds)

    ds = transform.coord_transform(ds)

    ds = aqdutils.magvar_correct(ds)

    ds = combine_vars(ds)

    if ds.attrs["VECBurstInterval"] != "CONTINUOUS":
        ds = dist_to_boundary(ds)

    ds = scale_analoginput(ds)

    # Shape into burst if not continuous mode
    if ds.attrs["VECBurstInterval"] != "CONTINUOUS":
        ds = reshape(ds)

    # Drop unused variables
    ds = ds_drop(ds)

    ds = qaqc.drop_vars(ds)

    # Rename DataArrays for EPIC compliance
    ds = aqdutils.ds_rename(ds)

    # Add EPIC and CMG attributes
    ds = aqdutils.ds_add_attrs(ds, inst_type="VEC")

    ds = associate_z_coord(ds)

    for var in ds.data_vars:
        # need to do this or else a "coordinates" attribute with value of "burst" hangs around
        ds[var].encoding["coordinates"] = None
        ds = qaqc.trim_min(ds, var)
        ds = qaqc.trim_max(ds, var)
        ds = qaqc.trim_min_diff(ds, var)
        ds = qaqc.trim_min_diff_pct(ds, var)
        ds = qaqc.trim_max_diff(ds, var)
        ds = qaqc.trim_max_diff_pct(ds, var)
        ds = qaqc.trim_med_diff(ds, var)
        ds = qaqc.trim_med_diff_pct(ds, var)
        ds = qaqc.trim_max_blip(ds, var)
        ds = qaqc.trim_max_blip_pct(ds, var)
        ds = qaqc.trim_bad_ens(ds, var)
        ds = qaqc.trim_bad_ens_indiv(ds, var)
        ds = qaqc.trim_fliers(ds, var)
        ds = qaqc.trim_warmup(ds, var)

    # after check for masking vars by other vars
    for var in ds.data_vars:
        ds = qaqc.trim_mask(ds, var)

    # Add start_time and stop_time attrs
    ds = utils.add_start_stop_time(ds)

    # Add history showing file used
    ds = utils.add_history(ds)

    ds = utils.add_standard_names(ds)

    ds["AGC_1202"].attrs.pop("standard_name")

    ds = ds.chunk({"time": 200000})

    ds = time_encoding(ds)

    if ds.attrs["VECBurstInterval"] != "CONTINUOUS":
        create_burst_nc(ds)

    if ds.attrs["VECBurstInterval"] == "CONTINUOUS":
        create_cont_nc(ds)

    return ds


def set_orientation(VEL):
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

    # User orientation refers to probe orientation
    # Nortek status code orientation refers to z-axis positive direction
    # See Nortek "The Comprehensive Manual - Velocimeters"
    # section 3.1.7 Orientation of Vector probes
    userorient = VEL.attrs["orientation"]
    # last bit of statuscode is orientation, which is saved as variable for vel transformations (inspried by dolfyn)
    sc = str(VEL["orientation"].isel(time=int(len(VEL["time"]) / 2)).values)
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

        warnings.warn(
            "Incorrect vector orientation will cause erroneous velocity transformations. Check deployment orientation details and vector manual to verify vector was oriented correctly"
        )

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
        "SNR1",
        "SNR2",
        "SNR3",
        "COR1",
        "COR2",
        "COR3",
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


def reshape(ds):

    t = ds["time"][ds["sample"].values == 1]

    for i in np.arange(0, len(t)):
        t2, samp = np.meshgrid(
            t[i],
            ds["sample"].sel(
                time=slice(
                    t[i],
                    t[i]
                    + np.timedelta64(ds.attrs["VECBurstInterval"], "s")
                    - np.timedelta64(
                        int(1 / ds.attrs["VECSamplingRate"] * 1000 * 0.5), "ms"
                    ),
                )
            ),
        )

        if i == 0:
            s = np.array(samp.transpose().flatten())
            t3 = np.array(t2.transpose().flatten())
        else:
            s = np.append(s, samp.transpose().flatten())
            t3 = np.append(t3, t2.transpose().flatten())

    ind = pd.MultiIndex.from_arrays((t3, s), names=("new_time", "new_sample"))

    ds = ds.sel(time=slice(t[0], ds["time"][-1])).assign(time=ind).unstack("time")

    ds = ds.drop("sample").rename({"new_time": "time", "new_sample": "sample"})

    return ds


def combine_vars(ds):

    ds["AGC_1202"] = (ds["AMP1"] + ds["AMP2"] + ds["AMP3"]) / 3

    veldim = ds.attrs["VECCoordinateSystem"].lower()

    if veldim == "xyz":
        veldim = "inst"

    ds["vel"] = xr.DataArray([ds.VEL1, ds.VEL2, ds.VEL3], dims=[veldim, "time"])
    ds["vel"].attrs.update(
        {
            "units": "m s-1",
            "long_name": f"Velocity, {veldim} coordinate system",
        }
    )

    ds["cor"] = xr.DataArray([ds.COR1, ds.COR2, ds.COR3], dims=["beam", "time"]).astype(
        "int32"
    )
    ds["amp"] = xr.DataArray([ds.AMP1, ds.AMP2, ds.AMP3], dims=["beam", "time"]).astype(
        "int32"
    )
    ds["snr"] = xr.DataArray([ds.SNR1, ds.SNR2, ds.SNR3], dims=["beam", "time"])

    return ds


def time_encoding(ds):
    """ensure we don't set dtypes uint for CF compliance"""
    if "units" in ds["time"].encoding:
        ds["time"].encoding.pop("units")

    # use time step to select time encoding
    tstep = ds["time"][1] - ds["time"][0]

    if tstep < np.timedelta64(1, "m"):

        histtext = f"make time encoding to dtype double because tstep {tstep} seconds is < 1 minute, round to milliseconds first"
        ds = utils.insert_history(ds, histtext)

        # round time to milliseconds first
        ds["time"] = ds["time"].dt.round("ms")
        ds["time"].encoding["dtype"] = "double"

    else:
        histtext = f"make time encoding int because tstep {tstep} seconds is >= 1 minute, round time to seconds first"
        ds = utils.insert_history(ds, histtext)

        # round time to seconds if time interval >= 1 minute
        ds["time"] = ds["time"].dt.round("s")

        # check time to make sure it fits in int32, assume seconds for time units
        utils.check_time_fits_in_int32(ds, "time")

        ds["time"].encoding["dtype"] = "i4"

    return ds


def create_burst_nc(ds):

    if "prefix" in ds.attrs:
        nc_filename = ds.attrs["prefix"] + ds.attrs["filename"] + "b-cal.nc"
    else:
        nc_filename = ds.attrs["filename"] + "b-cal.nc"

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
        nc_filename = dsmean.attrs["prefix"] + dsmean.attrs["filename"] + "-a.nc"
    else:
        nc_filename = dsmean.attrs["filename"] + "-a.nc"

    dsmean.to_netcdf(nc_filename)
    utils.check_compliance(nc_filename, conventions=dsmean.attrs["Conventions"])

    print("Done writing netCDF file", nc_filename)


def create_cont_nc(ds):

    if "prefix" in ds.attrs:
        nc_filename = ds.attrs["prefix"] + ds.attrs["filename"] + "cont-cal.nc"
    else:
        nc_filename = ds.attrs["filename"] + "cont-cal.nc"

    delayed_obj = ds.to_netcdf(nc_filename, compute=False)
    with ProgressBar():
        delayed_obj.compute()
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])

    print("Done writing netCDF file", nc_filename)
