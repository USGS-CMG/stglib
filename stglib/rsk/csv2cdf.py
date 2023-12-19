import warnings

import numpy as np
import pandas as pd
import xarray as xr
import yaml

from ..core import utils


def csv_to_cdf(metadata):
    basefile = metadata["basefile"]

    with open(basefile + "_metadata.txt") as f:
        meta = yaml.safe_load(f)

    print(f"Reading {basefile}_data.txt")
    try:
        df = pd.read_csv(basefile + "_data.txt", engine="pyarrow")
    except ValueError:
        warnings.warn(
            "Failed to load file using pyarrow; falling back to the slower default loader. Are you sure pyarrow is installed?"
        )
        df = pd.read_csv(
            basefile + "_data.txt"
        )  # pyarrow seems to fail on Python 3.7 on GitHub Actions

    df = df.rename(columns={"Time": "time"}).set_index("time")

    ds = df.to_xarray()

    # write out metadata first, then deal exclusively with xarray attrs
    ds = utils.write_metadata(ds, metadata)

    del metadata

    ds = utils.ensure_cf(ds)

    ds["time"] = pd.DatetimeIndex(ds["time"])

    ds = replace_spaces_in_var_names(ds)

    ds = rename_vars(ds, meta)

    ds = drop_unused_vars(ds)

    ds = set_up_instrument_and_sampling_attrs(ds, meta)

    is_profile = (
        (ds.attrs["sample_mode"] == "CONTINUOUS")
        and ("featureType" in ds.attrs)
        and (ds.attrs["featureType"] == "profile")
    )

    if ds.attrs["sample_mode"] == "WAVE":
        burst = pd.read_csv(basefile + "_burst.txt")

        burst = burst.rename(columns={"Time": "time"}).set_index("time")

        burst = burst.to_xarray()

        r = np.shape(burst.Pressure)[0]
        mod = r % ds.attrs["samples_per_burst"]
        if mod:
            print(
                "Number of rows is not a multiple of samples_per_burst; truncating to last full burst"
            )
            burst = burst.sel(time=burst.time[0:-mod])

        dsburst = xr.Dataset()

        dsburst["time"] = xr.DataArray(
            burst.time[0 :: int(ds.attrs["samples_per_burst"])].values, dims="time"
        )
        dsburst["time"] = pd.DatetimeIndex(dsburst["time"])

        dsburst["sample"] = xr.DataArray(
            range(ds.attrs["samples_per_burst"]), dims="sample"
        )

        for v in burst.data_vars:
            dsburst[v] = xr.DataArray(
                np.reshape(burst[v].values, (-1, int(ds.attrs["samples_per_burst"]))),
                dims=["time", "sample"],
            )

        if "Pressure" in dsburst:
            dsburst = get_metadata(dsburst, "Pressure", meta, "burstheader")
            dsburst = dsburst.rename({"Pressure": "P_1"})

        if "Burst" in dsburst:
            dsburst = dsburst.rename({"Burst": "burst"})
            dsburst["burst"] = dsburst["burst"].astype("int")

        if "Wave" in dsburst:
            dsburst = dsburst.drop("Wave")

        ds = ds.sel(time=dsburst.sel(sample=0).time)
        # sample gets added with the .sel above
        ds = ds.drop("sample")
        # drop the burst average value to replae with burst data
        ds = ds.drop("P_1")
        ds = xr.merge([ds, dsburst])

    elif is_profile:
        # work with profiles, e.g. CTD casts

        events = pd.read_csv(basefile + "_events.txt")
        events.rename(columns={"Time": "time"}, inplace=True)
        events["time"] = pd.to_datetime(events["time"])
        events.set_index("time", inplace=True)
        st = ["started" in x for x in events["Type"]]
        en = ["paused" in x for x in events["Type"]]
        starts = events[st].index
        ends = events[en].index
        # sometimes the first start event seems to be missing from the events file
        if starts[0] > ends[0]:
            starts = np.insert(starts, 0, ds.time[0].values)
        if starts.shape != ends.shape:
            raise ValueError("starts shape does not equal ends shape")

        pr = xr.Dataset()

        pr["profile"] = xr.DataArray(range(len(starts)), dims="profile")
        pr["profile"].attrs["cf_role"] = "profile_id"
        pr["profile"].encoding["dtype"] = "i4"

        timeavg = []
        rowsize = []
        for s, e in zip(starts, ends):
            dss = ds.time.sel(time=slice(s, e))
            timeavg.append(dss[0].values)
            rowsize.append(len(dss))

        pr["time"] = xr.DataArray(timeavg, dims="profile")
        pr["rowSize"] = xr.DataArray(rowsize, dims="profile")
        pr["rowSize"].attrs["long_name"] = "number of obs for this profile"
        pr["rowSize"].attrs["sample_dimension"] = "obs"
        pr["rowSize"].encoding["dtype"] = "i4"

        if "latitude" in ds.attrs and "longitude" in ds.attrs:
            if len(ds.attrs["latitude"]) == len(rowsize) and len(
                ds.attrs["longitude"]
            ) == len(rowsize):
                ds["latitude"] = xr.DataArray(
                    np.array(ds.attrs["latitude"]).astype(float), dims="profile"
                )
                ds["latitude"].attrs.update(
                    {
                        "units": "degree_north",
                        "axis": "Y",
                        "standard_name": "latitude",
                    }
                )

                ds["longitude"] = xr.DataArray(
                    np.array(ds.attrs["longitude"]).astype(float), dims="profile"
                )
                ds["longitude"].attrs.update(
                    {
                        "units": "degree_east",
                        "axis": "X",
                        "standard_name": "longitude",
                    }
                )

                ds.attrs.pop("latitude")
                ds.attrs.pop("longitude")
            else:
                raise ValueError(
                    f"size of latitude ({len(ds.attrs['latitude'])}) and longitude ({len(ds.attrs['longitude'])}) does not match number of profiles ({len(rowsize)})"
                )

        # dscp = ds.copy(deep=True)
        ds["obs"] = xr.DataArray(range(len(ds["time"])), dims="obs")
        ds["obs"].encoding["dtype"] = "i4"

        ds = ds.drop("time")

        ds = ds.rename({"obs": "time"}).set_coords("time").rename({"time": "obs"})

        ds = xr.merge([ds, pr])

    """
    # Set burst interval, [sec], USER DEFINED in instrument attr for continuous mode sampling
    if (ds.attrs["sample_mode"] == "CONTINUOUS") and ("wave_interval" in ds.attrs):
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

        ds["timenew"] = xr.DataArray(
            ds.time[0 :: int(ds.attrs["samples_per_burst"])].values, dims="timenew"
        )

        ds["sample"] = xr.DataArray(range(ds.attrs["samples_per_burst"]), dims="sample")

        for v in ["P_1", "T_28"]:
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
    """

    ds = utils.shift_time(ds, 0)

    # configure file
    cdf_filename = ds.attrs["filename"] + "-raw.cdf"

    if is_profile:
        ds.to_netcdf(cdf_filename, unlimited_dims=["obs"])
    else:
        ds.to_netcdf(cdf_filename, unlimited_dims=["time"])

    print("Finished writing data to %s" % cdf_filename)

    return ds


def create_lat_lon_vars_from_attrs(ds):
    ds["latitude"] = xr.DataArray(
        [ds.attrs["latitude"]],
        dims="latitude",
        attrs={"units": "degree_north", "standard_name": "latitude", "axis": "Y"},
    )

    ds["longitude"] = xr.DataArray(
        [ds.attrs["longitude"]],
        dims="longitude",
        attrs={"units": "degree_east", "standard_name": "longitude", "axis": "X"},
    )

    return ds


def replace_spaces_in_var_names(ds):
    for k in ds.data_vars:
        if " " in k:
            ds = ds.rename({k: k.replace(" ", "_")})

    return ds


def rename_vars(ds, meta):
    if "Turbidity" in ds:
        ds = get_metadata(ds, "Turbidity", meta)
        ds = ds.rename({"Turbidity": "Turb"})

        # ds["Turb"].attrs.update(
        #     {"units": "Nephelometric turbidity units (NTU)", "long_name": "Turbidity"}
        # )

    if "Pressure" in ds:
        ds = get_metadata(ds, "Pressure", meta)
        ds = ds.rename({"Pressure": "P_1"})

    if "Depth" in ds:
        ds = ds.drop("Depth")

    if "Conductivity" in ds:
        ds = get_metadata(ds, "Conductivity", meta)
        ds = ds.rename({"Conductivity": "C_51"})

    if "Salinity" in ds:
        ds = get_metadata(ds, "Salinity", meta)
        ds = ds.rename({"Salinity": "S_41"})

    if "Temperature" in ds:
        ds = get_metadata(ds, "Temperature", meta)
        ds = ds.rename({"Temperature": "T_28"})

    if "Specific_conductivity" in ds:
        ds = get_metadata(ds, "Specific_conductivity", meta)
        ds = ds.rename({"Specific_conductivity": "SpC_48"})

    return ds


def drop_unused_vars(ds):
    for var in [
        "Sea_pressure",
        "Depth",
        "Speed_of_sound",
        "Density_anomaly",
        "Tidal_slope",
    ]:
        if var in ds:
            ds = ds.drop(var)

    return ds


def set_up_instrument_and_sampling_attrs(ds, meta):
    ds.attrs["serial_number"] = str(meta["instrument"]["serial"])
    ds.attrs["instrument_type"] = str(meta["instrument"]["model"])

    for v in ["fwtype", "fwversion"]:
        ds.attrs[v] = str(meta["instrument"][v])

    ds.attrs["sample_mode"] = meta["sampling"]["mode"]

    ds.attrs["sample_interval"] = int(meta["sampling"]["period"]) / 1000

    if "burstinterval" in meta["sampling"]:
        ds.attrs["burst_interval"] = meta["sampling"]["burstinterval"] / 1000

    if "burstcount" in meta["sampling"]:
        ds.attrs["samples_per_burst"] = meta["sampling"]["burstcount"]

    return ds


def get_metadata(ds, var, meta, field="dataheader"):
    dh = meta[field]
    metanames = np.array([x["name"] for x in dh])
    thisname = np.where(metanames == var.replace("_", " "))[0][0]
    for v in ["units", "calibration", "ranging"]:
        if v in dh[thisname]:
            if type(dh[thisname][v]) == dict:
                ds[var].attrs[v] = str(dh[thisname][v])
            else:
                ds[var].attrs[v] = dh[thisname][v]

    return ds
