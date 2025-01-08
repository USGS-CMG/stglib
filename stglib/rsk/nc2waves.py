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

    for k in ["wp_peak", "wh_4061", "wp_4060", "pspec"]:
        ds[k] = spec[k]

    # ds = utils.create_water_depth(ds)

    ds = utils.create_water_depth_var(ds)
    ds = utils.create_water_level_var(ds)

    for k in ["P_1", "P_1ac", "sample", "T_28", "water_level_filt"]:
        if k in ds:
            ds = ds.drop_vars(k)

    ds = qaqc.drop_vars(ds)

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


def nc_to_diwasp(nc_filename):
    """
    Process burst data to make DIWASP derived wave statistics and spectra
    """
    print("starting RBR nc2diwasp")
    ds = xr.open_dataset(nc_filename)

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
        ds = make_wave_bursts(ds)

    if "sample_rate" not in ds.attrs:
        ds.attrs["sample_rate"] = 1 / ds.attrs["sample_interval"]

    if "diwasp" not in ds.attrs:
        # 'pres' is only option for a single pressure logger
        ds.attrs["diwasp"] = "pres"

    if "diwasp" in ds.attrs:
        data_type = ds.attrs["diwasp"]
        if data_type not in ["pres"]:
            raise ValueError(
                f"data type {data_type} is not recognized for RBR pressure logger current options are ['pres']"
            )

        print(f"Running DIWASP using {ds.attrs['diwasp']} input data")

        pz = ds.attrs["initial_instrument_height"]
        layout = np.atleast_2d([0, 0, pz]).T
        diwasp = waves.make_diwasp_elev_pres(ds, data_type=data_type, layout=layout)
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
    ]:

        if k in diwasp:
            ds[k] = diwasp[k]

    # add diwasp attrs
    for k in diwasp.attrs:
        ds.attrs[k] = diwasp.attrs[k]

    # rename Fspec based on input datatype
    ds = utils.rename_diwasp_fspec(ds)

    ds = utils.create_water_depth_var(ds)

    for k in [
        "burst",
        "sample",
        "P_1",
        "P_1ac",
        "Tx_1211",
    ]:
        if k in ds:
            ds = ds.drop_vars(k)

    # ds = qaqc.drop_vars(ds)

    # use "epic" names
    ds.attrs["diwasp_names"] = "epic"
    ds = utils.rename_diwasp_wave_vars(ds)

    ds = drop_unused_dims(ds)

    ds = utils.trim_max_wp(ds)

    ds = utils.trim_min_wh(ds)

    ds = utils.trim_max_wh(ds)

    ds = utils.trim_wp_ratio(ds)

    # Add attrs
    ds = utils.ds_add_wave_attrs(ds)

    # assign min/max (need to do this after trimming):
    ds = utils.add_min_max(ds)

    # round time to minutes to make sure fits in dtype i4. will be fine for wave burst start times
    ds["time"] = ds["time"].dt.round("min")

    nc_filename = ds.attrs["filename"] + "s_diwasp-a.nc"

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
