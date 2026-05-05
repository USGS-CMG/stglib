import re
import time
import warnings

import numpy as np
import xarray as xr

from ..core import qaqc, utils, waves
from ..sig.cdf2nc import add_water_level
from ..sig.nc2waves import (
    drop_unused_dims,
    make_wave_bursts_mi,
    make_wave_vars,
    make_waves_vdims,
)


def nc_to_waves(nc_filename, salwtemp=None):
    """
    Process burst data to wave statistics
    """

    ds = xr.load_dataset(nc_filename)

    spec = waves.make_waves_ds(ds)

    for k in ["wp_peak", "wh_4061", "wp_4060", "pspec"]:
        ds[k] = spec[k]

    dopuv = False
    if "puv" in ds.attrs:
        if ds.attrs["puv"].lower() == "true":
            dopuv = True

    if dopuv:
        ds = do_puv(ds)

    # keep burst mean P_1 and P_1ac for reference
    for k in ["P_1", "P_1ac"]:
        if k in ds:
            ds[k] = ds[k].mean(dim="sample", keep_attrs=True)

    if "pressure_sensor_height" in ds.attrs:
        ds = utils.create_water_depth_var(
            ds.mean(dim="sample", keep_attrs=True),
            psh="pressure_sensor_height",
            salwtemp=salwtemp,
        )
    else:
        ds = utils.create_water_depth_var(
            ds.mean(dim="sample", keep_attrs=True), salwtemp=salwtemp
        )

    ds = add_water_level(ds, salwtemp=salwtemp)

    for k in [
        "P_1",
        "P_1ac",
        "burst",
        "brange",
        "sample",
        "u_1205",
        "v_1206",
        "w_1204",
        "vel1_1277",
        "vel2_1278",
        "vel3_1279",
        "vel",
        "vrange",
        "amp",
        "amp_avg",
        "AGC_1202",
        "AGC1_1221",
        "AGC2_1222",
        "AGC3_1223",
        "SNR1",
        "SNR2",
        "SNR3",
        "snr",
        "cor1_1285",
        "cor2_1286",
        "cor3_1287",
        "cor",
        "AnalogInput1",
        "AnalogInput2",
        "Hdg_1215",
        "Ptch_1216",
        "Roll_1217",
        "Bat_106",
        "Tx_1211",
        "SV_80",
        "orientation",
        "orientmat",
        "TransMatrix",
    ]:
        if k in ds:
            ds = ds.drop_vars(k)

    ds = make_waves_vdims(ds, wtype="nc2waves")

    ds = drop_unused_dims(ds)

    ds = qaqc.drop_vars(ds)

    ds = utils.trim_max_wp(ds)

    ds = utils.trim_min_wh(ds)

    ds = utils.trim_max_wh(ds)

    ds = utils.trim_wp_ratio(ds)

    if dopuv:
        ds = waves.puv_qaqc(ds)

    # Add attrs
    ds = utils.ds_add_wave_attrs(ds)

    # assign min/max (need to do this after trimming):
    ds = utils.add_min_max(ds)

    nc_filename = ds.attrs["filename"] + "s-a.nc"

    ds.to_netcdf(nc_filename, unlimited_dims=["time"])
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])
    print("Done writing netCDF file", nc_filename)

    return ds


def do_puv(ds):
    print("Running puv_quick")

    for k in ["pressure_sensor_height", "velocity_sample_volume_height"]:
        if k not in ds.attrs:
            raise KeyError(f"{k} must be specified to run PUV")

    N, M = np.shape(ds["u_1205"].squeeze())

    if "puv_first_frequency_cutoff" in ds.attrs:
        first_frequency_cutoff = ds.attrs["puv_first_frequency_cutoff"]
    else:
        first_frequency_cutoff = 1 / 10

    if "puv_last_frequency_cutoff" in ds.attrs:
        last_frequency_cutoff = ds.attrs["puv_last_frequency_cutoff"]
    else:
        last_frequency_cutoff = 1 / 2.5

    if "P_1ac" in ds:
        pvar = "P_1ac"
    else:
        pvar = "P_1"

    u = ds["u_1205"].squeeze().transpose("time", "sample")
    v = ds["v_1206"].squeeze().transpose("time", "sample")

    psh = ds.attrs["pressure_sensor_height"]
    huv = ds.attrs["velocity_sample_volume_height"]

    ds = waves.call_puv_quick_vectorized(
        ds,
        pvar=pvar,
        u=u,
        v=v,
        psh=psh,
        huv=huv,
        first_frequency_cutoff=first_frequency_cutoff,
        last_frequency_cutoff=last_frequency_cutoff,
    )

    return ds


def nc_to_diwasp(nc_filename, salwtemp=None):
    """
    Process burst data to make DIWASP derived wave statistics and spectra
    """

    start_time = time.time()
    ds = xr.open_dataset(nc_filename, chunks={"time": 200000})

    # check to see if need to make wave burst from continuous data
    if (ds.attrs["sample_mode"] == "CONTINUOUS") and ("wave_interval" in ds.attrs):
        # check for wave_start_time attrs
        if "wave_start_time" in ds.attrs:
            print(
                f"trimming continuous data to start at users specified wave start time {ds.attrs['wave_start_time']}"
            )
            ds = ds.sel(
                time=slice(np.datetime64(ds.attrs["wave_start_time"]), ds["time"][-1])
            )
        # make wave burst ncfile from continuous data if wave_interval is specified
        ds = make_wave_bursts_mi(ds)

    if "BURST" in ds.attrs["sample_mode"]:
        if "wave_interval" not in ds.attrs:
            if "burst_interval" in ds.attrs:
                ds.attrs["wave_interval"] = ds.attrs["burst_interval"]
            else:
                ds.attrs["wave_interval"] = int(
                    (ds["time"][1] - ds["time"][0]).dt.round("s").values
                    / np.timedelta64(1, "s")
                )
        # make data set for wave analysis
        ds = make_wave_vars(ds)

    # check to see if wave duration is specified and if so trim burst samples accordingly
    if "wave_duration" in ds.attrs:
        if ds.attrs["wave_duration"] <= ds.attrs["wave_interval"]:
            nsamps = int(ds.attrs["wave_duration"] * ds.attrs["sample_rate"])
            ds = ds.isel(sample=slice(0, nsamps))
            ds.attrs["wave_samples_per_burst"] = nsamps

        else:
            warnings.warn(
                f"wave_duration {ds.attrs.wave_duration} cannot be greater than wave_interval {ds.attrs.wave_interval}, this configuration setting will be ignored"
            )

    if "diwasp" not in ds.attrs:
        # if not specified choose 'suv' method for signature
        ds.attrs["diwasp"] = "puv"
        print("diwasp processing type was not specified, defaulting to 'puv' method")

    if "diwasp" in ds.attrs:
        data_type = ds.attrs["diwasp"]
        if data_type not in ["suv", "puv", "elev", "pres", "optimized", "optimized-nd"]:
            raise ValueError(
                f"data type {data_type} is not recognized current options are ['suv', 'puv', 'elev', 'pres', 'optimized', 'optimized-nd']"
            )

    # if "diwasp_ibin" in ds.attrs:
    #    ibin = ds.attrs["diwasp_ibin"]
    # else:
    ibin = 0

    if data_type == "puv":
        print(f"Running DIWASP using {data_type} input data")
        layout = make_diwasp_layout(ds, data_type=data_type, ibin=ibin)
        diwasp = waves.make_diwasp_ds(ds, data_type=data_type, layout=layout, ibin=ibin)
        ds = utils.ds_add_pydiwasp_history(ds)

    elif data_type == "pres":
        print(f"Running DIWASP using {ds.attrs['diwasp']} input data")
        layout = make_diwasp_layout(ds, data_type=data_type, ibin=ibin)
        diwasp = waves.make_diwasp_ds(ds, data_type=data_type, layout=layout, ibin=ibin)
        ds = utils.ds_add_pydiwasp_history(ds)

    else:
        raise ValueError(
            f"DIWASP input type {ds.attrs['diwasp']} is not currently supported for {ds.attrs['instument_type']} in stglib"
        )

    for k in [
        "diwasp_frequency",
        "diwasp_direction",
        "diwasp_hs",
        "diwasp_tp",
        "diwasp_tm",
        "diwasp_dtp",
        "diwasp_dp",
        "diwasp_dm",
        "diwasp_dspec",
        "diwasp_fspec",
        "diwasp_type",
    ]:

        if k in diwasp:
            ds[k] = diwasp[k]

    # add diwasp attrs
    for k in diwasp.attrs:
        ds.attrs[k] = diwasp.attrs[k]

    ds = utils.rename_diwasp_fspec(ds)

    if "pressure_sensor_height" in ds.attrs:
        ds = utils.create_water_depth_var(
            ds.mean(dim="sample", keep_attrs=True),
            psh="pressure_sensor_height",
            salwtemp=salwtemp,
        )
    else:
        ds = utils.create_water_depth_var(
            ds.mean(dim="sample", keep_attrs=True), salwtemp=salwtemp
        )

    ds = add_water_level(ds, salwtemp=salwtemp)

    for k in [
        "P_1",
        "P_1ac",
        "burst",
        "brange",
        "sample",
        "u_1205",
        "v_1206",
        "w_1204",
        "vel1_1277",
        "vel2_1278",
        "vel3_1279",
        "vel",
        "vrange",
        "amp",
        "amp_avg",
        "AGC_1202",
        "AGC1_1221",
        "AGC2_1222",
        "AGC3_1223",
        "SNR1",
        "SNR2",
        "SNR3",
        "snr",
        "cor1_1285",
        "cor2_1286",
        "cor3_1287",
        "cor",
        "AnalogInput1",
        "AnalogInput2",
        "Hdg_1215",
        "Ptch_1216",
        "Roll_1217",
        "Bat_106",
        "Tx_1211",
        "SV_80",
        "orientation",
        "orientmat",
        "TransMatrix",
    ]:

        if k in ds:
            ds = ds.drop_vars(k)

    ds = qaqc.drop_vars(ds)

    # use "epic" names
    ds.attrs["diwasp_names"] = "epic"
    ds = utils.rename_diwasp_wave_vars(ds)

    ds = make_waves_vdims(ds, wtype="nc2diwasp")

    ds = drop_unused_dims(ds)

    ds = utils.trim_max_wp(ds)

    ds = utils.trim_min_wh(ds)

    ds = utils.trim_max_wh(ds)

    ds = utils.trim_wp_ratio(ds)

    # Add attrs
    ds = utils.ds_add_wave_attrs(ds)

    # clean-up attrs
    ds = drop_attrs(ds)

    # assign min/max (need to do this after trimming):
    ds = utils.add_min_max(ds)

    ds = utils.add_delta_t(ds)

    # round time to minutes to make sure fits in dtype i4. will be fine for wave burst start times
    ds["time"] = ds["time"].dt.round("min")

    nc_filename = ds.attrs["filename"] + "s_diwasp-a.nc"

    ds.to_netcdf(nc_filename, unlimited_dims=["time"])
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])

    print("Done writing netCDF file", nc_filename)

    end_time = time.time()
    print(
        f"Processing nc2diwasp for {ds.attrs['instrument_type']} data  in deployment {ds.attrs['filename']} completed"
    )
    print(f"elapsed time = {end_time-start_time:.0f} seconds")

    return ds


def make_diwasp_layout(ds, data_type=None, ibin=None):
    """
    Make layout for DIWASP wave processing input
    """

    # check for pressure sensor x,y offsets in yaml file
    if "px_offset" in ds.attrs:
        px = ds.attrs["px_offset"]
    else:
        px = 0

    if "py_offset" in ds.attrs:
        py = ds.attrs["py_offset"]
    else:
        py = 0

    if data_type == "puv":
        # datatypes = ["pres", "velx", "vely"]
        sxyz = [px, py, ds.attrs["pressure_sensor_height"]]
        if ds.attrs["orientation"].lower() == "up":
            velz = ds.attrs["initial_instrument_height"] + ds["bindist"][ibin].values

        elif ds.attrs["orientation"].lower() == "down":
            velz = ds.attrs["initial_instrument_height"] - ds["bindist"][ibin].values

        uxyz = vxyz = [0, 0, velz]

        layout = np.array([sxyz, uxyz, vxyz]).T

    elif data_type == "pres":
        # datatypes=['elev']
        sxyz = np.atleast_2d([px, py, ds.attrs["pressure_sensor_height"]]).T
        layout = sxyz

    return layout


def drop_attrs(ds):
    """Drop some global attrs"""

    rm = []  # initialize list of attrs to be removed
    exclude = []  # initialize attrs to exclude

    for j in ds.attrs:
        if re.match("^VEC", j):
            rm.append(j)

    for k in rm:
        if k not in exclude:
            del ds.attrs[k]

    return ds
