import numpy as np
import pandas as pd
import xarray as xr

from ..core import qaqc, utils, waves


def nc_to_waves(nc_filename):
    """
    Process burst data to wave statistics
    """

    ds = xr.open_dataset(nc_filename)

    # Convert Hertz to sample interval in seconds
    ds.attrs["sample_interval"] = 1 / float(ds.attrs["WaveSampleRate"])

    # Check to see if need to make smaller wave bursts from really long wave bursts
    if "calculated_wave_interval" in ds.attrs:
        # Divide large burst into smaller bursts at specified calculated_wave_interval
        ds = make_wave_bursts(ds)

    spec = waves.make_waves_ds(ds)

    for k in ["wp_peak", "wh_4061", "wp_4060", "pspec"]:
        ds[k] = spec[k]

    ds = utils.create_water_depth_var(ds)

    # Drop unneeded variables/attributes
    for k in ["P_1", "P_1ac", "sample", "T_28"]:
        if k in ds:
            ds = ds.drop_vars(k)
    del ds.attrs["sample_interval"]  # Delete because info already in gatts

    ds = qaqc.drop_vars(ds)

    ds = utils.trim_max_wp(ds)

    ds = utils.trim_min_wh(ds)

    ds = utils.trim_max_wh(ds)

    ds = utils.trim_wp_ratio(ds)

    # Add attrs
    ds = utils.ds_add_wave_attrs(ds)

    # assign min/max (need to do this after trimming):
    ds = utils.add_min_max(ds)

    nc_filename = ds.attrs["filename"] + "w-cal.nc"

    ds.to_netcdf(nc_filename, unlimited_dims=["time"])
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])
    print("Done writing netCDF file", nc_filename)

    return ds


def make_wave_bursts(ds):
    ds.attrs["samples_per_burst"] = int(
        ds.attrs["calculated_wave_interval"] / ds.attrs["sample_interval"]
    )

    # Multiply dimensions of P_1 to get total number of samples
    # no_of_samples = np.shape(ds.P_1)[0] * np.shape(ds.P_1)[1]

    # Calculate remaining rows if number of samples is not a multiple of bursts
    # rem = no_of_samples % ds.attrs["samples_per_burst"]
    # if rem:
    #     print(
    #         "Number of rows is not a multiple of samples_per_burst; truncating to last full burst"
    #     )
    #     ds = ds.sel(time=ds.time[0:-rem])

    ds.time.encoding.pop("dtype")

    # Calculate how many rows to subdivide each burst
    rows = float(ds.attrs["BurstDuration"]) / float(
        ds.attrs["calculated_wave_interval"]
    )

    # Calculate time difference for each row
    delta_t = int(float(ds.attrs["TideInterval"]) / rows)

    last_time = ds.time[-1].values + np.timedelta64(
        int(float(ds.attrs["TideInterval"]) - delta_t), "m"
    )

    delta_t = str(delta_t) + "min"

    new_date_rng = pd.date_range(start=ds.time[0].values, end=last_time, freq=delta_t)

    ds["new_time"] = xr.DataArray(new_date_rng, dims="new_time")

    ds["new_sample"] = xr.DataArray(
        range(ds.attrs["samples_per_burst"]), dims="new_sample"
    )

    for v in ["P_1", "P_1ac"]:
        if v in ds:
            attrsbak = ds[v].attrs
            ds[v] = xr.DataArray(
                np.reshape(ds[v].values, (-1, int(ds.attrs["samples_per_burst"]))),
                dims=["new_time", "new_sample"],
            )
            ds[v].attrs = attrsbak

    ds = ds.rename({"time": "timeold"})
    ds = ds.drop("timeold")

    ds = ds.rename({"sample": "sampleold"})
    ds = ds.drop("sampleold")

    ds = ds.rename({"new_sample": "sample"})
    ds = ds.rename({"new_time": "time"})

    return ds
