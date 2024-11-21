import re

import numpy as np
import xarray as xr

from ..core import qaqc, utils, waves


def nc_to_waves(nc_filename):
    """
    Process burst data to wave statistics
    """

    ds = xr.open_dataset(nc_filename)

    # check to see if need to make wave burst from continuous data
    if (ds.attrs["sample_mode"] == "CONTINUOUS") and ("wave_interval" in ds.attrs):
        # make wave burst ncfile from continuous data if wave_interval is specified
        ds = make_wave_bursts(ds)

    spec = waves.make_waves_ds(ds)
    ds = utils.ds_add_waves_history(ds)

    for k in ["wp_peak", "wh_4061", "wp_4060", "pspec"]:
        ds[k] = spec[k]

    dodiwasp = False
    if "diwasp" in ds.attrs:
        print("Running DIWASP")
        dodiwasp = True
        diwasp = waves.make_diwasp_puv_suv(ds, freqs=spec["frequency"].values)
        ds = utils.ds_add_pydiwasp_history(ds)

    for k in ["Hs", "Tp", "DTp", "Dp", "Dm", "dspec"]:
        ds[f"diwasp_{k}"] = diwasp[k]

    # add diwasp attrs
    for k in diwasp.attrs:
        ds.attrs[k] = diwasp.attrs[k]

    # ds = utils.create_water_depth(ds)

    ds = utils.create_water_depth_var(ds)

    for k in [
        "burst",
        "sample",
        "P_1",
        "P_1ac",
        "u_1205",
        "v_1206",
        "w_1204",
        "Tx_1211",
        "SV_80",
        "vel",
        "amp",
        "cor",
        "brangeAST",
        "TransMatrix",
    ]:
        if k in ds:
            ds = ds.drop_vars(k)

    # ds = qaqc.drop_vars(ds)

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

    # round time to minutes to make sure fits in dtype i4. will be fine for wave burst start times
    ds["time"] = ds["time"].dt.round("min")

    nc_filename = ds.attrs["filename"] + "s-a.nc"

    ds.to_netcdf(nc_filename, unlimited_dims=["time"])
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])

    """
    if dodiwasp:
        diwasp_filename = ds.attrs["filename"] + "_diwasp-cal.nc"
        diwasp.to_netcdf(diwasp_filename)
        utils.check_compliance(diwasp_filename, conventions=ds.attrs["Conventions"])
    """

    print("Done writing netCDF file", nc_filename)

    return ds


def make_wave_bursts(ds):
    # wave_interval is [sec] interval for wave statistics for continuous data
    ds.attrs["wave_samples_per_burst"] = int(
        ds.attrs["wave_interval"] / ds.attrs["sample_interval"]
    )
    # burst_interval is equivalent to wave_interval [sec]
    # ds.attrs["wave_burst_interval"] = ds.attrs["wave_interval"]
    # burst_length is the number of data points in the burst
    # ds.attrs["wave_burst_length"] = ds.attrs["wave_samples_per_burst"]
    r = np.shape(ds.P_1)[0]
    mod = r % ds.attrs["wave_samples_per_burst"]
    if mod:
        print(
            "Number of rows is not a multiple of samples_per_burst; truncating to last full burst"
        )
        ds = ds.sel(time=ds.time[0:-mod])

    ds.time.encoding.pop("dtype")

    ds["timenew"] = xr.DataArray(
        ds.time[0 :: int(ds.attrs["wave_samples_per_burst"])].values, dims="timenew"
    )

    ds["samplenew"] = xr.DataArray(
        range(ds.attrs["wave_samples_per_burst"]), dims="samplenew"
    )

    for v in ["P_1", "P_1ac", "Tx_1211", "brangeAST"]:
        if v in ds:
            attrsbak = ds[v].attrs
            ds[v] = xr.DataArray(
                np.reshape(ds[v].values, (-1, int(ds.attrs["wave_samples_per_burst"]))),
                dims=["timenew", "samplenew"],
            )
            ds[v].attrs = attrsbak

    for v in ["u_1205", "v_1206", "w_1204"]:
        if v in ds:
            attrsbak = ds[v].attrs
            # check dims
            if ds[v].dims == ("time", "z"):
                ds[v] = xr.DataArray(
                    np.reshape(
                        ds[v].values,
                        (-1, int(ds.attrs["wave_samples_per_burst"]), len(ds["z"])),
                    ).transpose(0, 2, 1),
                    dims=["timenew", "z", "samplenew"],
                )
                # u=np.reshape(ds['u_1205'].values, (-1,int(ds.attrs["samples_per_burst"]), len(ds['z']))).transpose(0,2,1)
                ds[v].attrs = attrsbak
            else:
                raise ValueError(
                    f"Not able to apply median filter because kernel size specified {kernel_size} is not an odd whole number"
                )

    for v in ["vel", "cor", "amp"]:
        if v in ds:
            attrsbak = ds[v].attrs
            ds[v] = xr.DataArray(
                np.reshape(
                    ds[v].values,
                    (
                        len(ds["beam"]),
                        -1,
                        int(ds.attrs["wave_samples_per_burst"]),
                        len(ds["z"]),
                    ),
                ).transpose(0, 1, 3, 2),
                dims=["beam", "timenew", "z", "samplenew"],
            )
            # u=np.reshape(ds['u_1205'].values, (-1,int(ds.attrs["samples_per_burst"]), len(ds['z']))).transpose(0,2,1)
            ds[v].attrs = attrsbak

    ds = ds.rename({"time": "timeold"})
    ds = ds.rename({"timenew": "time"})
    ds = ds.drop_dims("timeold")

    if "sample" in ds:
        ds = ds.rename({"sample": "sampleold"})
        ds = ds.rename({"samplenew": "sample"})
        ds = ds.drop_vars("sampleold")
    else:
        ds = ds.rename({"samplenew": "sample"})

    return ds


def drop_unused_dims(ds):
    """only keep dims that will be in the final files"""
    thedims = []
    for v in ds.data_vars:
        for x in ds[v].dims:
            thedims.append(x)

    for x in ds.dims:
        if x not in thedims:
            ds = ds.drop_vars(x)

    return ds


def drop_attrs(ds):
    """Drop some global attrs"""

    rm = []  # initialize list of attrs to be removed
    exclude = []  # initialize attrs to exclude

    for j in ds.attrs:
        if re.match(f"^SIG", j):
            rm.append(j)

    for k in rm:
        if k not in exclude:
            del ds.attrs[k]

    return ds
