import numpy as np
import xarray as xr

from ..core import utils, waves


def nc_to_waves(nc_filename):
    """
    Process burst data to wave statistics
    """

    ds = utils.open_time_2d_dataset(nc_filename)  # this will deal with a cf file, too

    # check to see if need to make wave burst from continuous data
    if (ds.attrs["sample_mode"] == "CONTINUOUS") and ("wave_interval" in ds.attrs):
        # make wave burst ncfile from continuous data if wave_interval is specified
        ds = make_wave_bursts(ds)

    spec = waves.make_waves_ds(ds)

    for k in ["wp_peak", "wh_4061", "wp_4060", "pspec"]:
        ds[k] = spec[k]

    # ds = utils.create_water_depth(ds)

    ds = utils.create_water_depth_var(ds)

    for k in ["P_1", "P_1ac", "sample", "T_28"]:
        if k in ds:
            ds = ds.drop_vars(k)

    ds = utils.trim_max_wp(ds)

    ds = utils.trim_min_wh(ds)

    ds = utils.trim_max_wh(ds)

    ds = utils.trim_wp_ratio(ds)

    # Add attrs
    ds = utils.ds_add_wave_attrs(ds)

    # assign min/max (need to do this after trimming):
    ds = utils.add_min_max(ds)

    nc_filename = ds.attrs["filename"] + "s-a.nc"

    ds.to_netcdf(nc_filename, unlimited_dims=["time"])
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])
    print("Done writing netCDF file", nc_filename)

    return ds


def make_wave_bursts(ds):
    # wave_interval is [sec] interval for wave statistics for continuous data
    ds.attrs["samples_per_burst"] = int(
        ds.attrs["wave_interval"] / ds.attrs["sample_interval"]
    )
    # burst_interval is equivalent to wave_interval [sec]
    ds.attrs["burst_interval"] = ds.attrs["wave_interval"]
    # burst_length is the number of data points in the burst
    ds.attrs["burst_length"] = ds.attrs["samples_per_burst"]
    r = np.shape(ds.P_1)[0]
    mod = r % ds.attrs["samples_per_burst"]
    if mod:
        print(
            "Number of rows is not a multiple of samples_per_burst; truncating to last full burst"
        )
        ds = ds.sel(time=ds.time[0:-mod])

    ds.time.encoding.pop("dtype")

    ds["timenew"] = xr.DataArray(
        ds.time[0 :: int(ds.attrs["samples_per_burst"])].values, dims="timenew"
    )

    ds["sample"] = xr.DataArray(range(ds.attrs["samples_per_burst"]), dims="sample")

    for v in ["P_1", "P_1ac", "T_28"]:
        if v in ds:
            attrsbak = ds[v].attrs
            ds[v] = xr.DataArray(
                np.reshape(ds[v].values, (-1, int(ds.attrs["samples_per_burst"]))),
                dims=["timenew", "sample"],
            )
            ds[v].attrs = attrsbak

    ds = ds.rename({"time": "timeold"})
    ds = ds.rename({"timenew": "time"})
    ds = ds.drop("timeold")

    return ds
