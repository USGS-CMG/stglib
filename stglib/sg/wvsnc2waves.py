import numpy as np
import pandas as pd
import xarray as xr

from ..core import qaqc, utils, waves


def nc_to_waves(nc_filename):
    """
    Process burst data to wave statistics
    """

    ds = xr.open_dataset(nc_filename)

    # Check to see if need to make smaller wave bursts from really long wave bursts
    if "calculated_wave_interval" in ds.attrs:
        # Divide large burst into smaller bursts at specified calculated_wave_interval
        ds = make_wave_bursts(ds)

    spec = waves.make_waves_ds(ds)

    for k in ["wp_peak", "wh_4061", "wp_4060", "pspec"]:
        ds[k] = spec[k]

    ds = utils.create_water_level_var(ds)
    ds = utils.create_water_depth_var(ds)

    # Drop unneeded variables/attributes
    for k in ["P_1", "P_1ac", "sample", "T_28"]:
        if k in ds:
            ds = ds.drop_vars(k)

    # Call all QAQC
    ds = qaqc.call_qaqc(ds)
    ds = utils.trim_max_wp(ds)
    ds = utils.trim_min_wh(ds)
    ds = utils.trim_max_wh(ds)
    ds = utils.trim_wp_ratio(ds)

    # Clip
    ds = utils.clip_ds(ds, wvs=True)

    # Add attrs
    ds = utils.ds_add_wave_attrs(ds)

    # assign min/max (need to do this after trimming):
    ds = utils.add_min_max(ds)
    ds = utils.add_delta_t(ds)

    nc_filename = ds.attrs["filename"] + "s-a.nc"

    ds.to_netcdf(nc_filename, unlimited_dims=["time"])
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])
    print("Done writing netCDF file", nc_filename)

    return ds


def make_wave_bursts(ds):
    ds.attrs["samples_per_burst"] = int(
        ds.attrs["calculated_wave_interval"] / ds.attrs["sample_interval"]
    )

    # Calculate how many rows to subdivide each burst
    rows = float(ds.attrs["SGBurstDuration"]) / float(
        ds.attrs["calculated_wave_interval"]
    )

    # Define time interval
    delta_t = int(ds.attrs["calculated_wave_interval"])
    delta_t = f"{delta_t}s"

    timestamp = []
    for t in range(0, ds["time"].size):

        # Make timestamp
        new_time = np.array(
            pd.date_range(start=ds.time[t].values, periods=rows, freq=delta_t)
        )
        timestamp.append(new_time)

    # Append all timestamps to array
    ds["date"] = np.concatenate(timestamp)

    # new_date_rng = pd.date_range(start=ds.time[0].values, end=last_time, freq=delta_t)

    # ds["new_time"] = xr.DataArray(new_date_rng, dims="new_time")

    ds["new_sample"] = xr.DataArray(
        np.arange(ds.attrs["samples_per_burst"], dtype=np.int32),
        dims="new_sample",
    )

    for v in ["P_1", "P_1ac"]:
        if v in ds:
            attrsbak = ds[v].attrs
            ds[v] = xr.DataArray(
                np.reshape(ds[v].values, (-1, int(ds.attrs["samples_per_burst"]))),
                dims=["date", "new_sample"],
            )
            ds[v].attrs = attrsbak

    ds = ds.rename({"time": "timeold"})
    ds = ds.drop("timeold")

    ds = ds.rename({"sample": "sampleold"})
    ds = ds.drop("sampleold")

    ds = ds.rename({"new_sample": "sample"})
    ds = ds.rename({"date": "time"})

    return ds
