import pandas as pd
import yaml
import numpy as np
import xarray as xr
from ..core import utils


def csv_to_cdf(metadata):

    basefile = metadata["basefile"]

    with open("_".join(basefile.split("_")[0:-1]) + "_metadata.txt") as f:
        meta = yaml.safe_load(f)

    df = pd.read_csv(basefile + ".txt", infer_datetime_format=True)

    df = df.rename(columns={"Time": "time"}).set_index("time")

    ds = df.to_xarray()

    # write out metadata first, then deal exclusively with xarray attrs
    ds = utils.write_metadata(ds, metadata)

    del metadata

    ds["time"] = pd.DatetimeIndex(ds["time"])

    ds.attrs["serial_number"] = meta["instrument"]["serial"]
    for v in ["model", "fwtype", "fwversion"]:
        ds.attrs[v] = meta["instrument"][v]
    ds.attrs["sample_mode"] = meta["sampling"]["mode"]

    for k in ds.data_vars:
        if " " in k:
            ds = ds.rename({k: k.replace(" ", "_")})

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

    for var in [
        "Sea_pressure",
        "Depth",
        "Speed_of_sound",
        "Density_anomaly",
        "Tidal_slope",
    ]:
        if var in ds:
            ds = ds.drop(var)

    ds.attrs["sample_interval"] = int(meta["sampling"]["period"]) / 1000
    if "burstinterval" in meta["sampling"]:
        ds.attrs["burst_interval"] = meta["sampling"]["burstinterval"] / 1000
    if "burstcount" in meta["sampling"]:
        ds.attrs["samples_per_burst"] = meta["sampling"]["burstcount"]

    if ds.attrs["sample_mode"] == "WAVE":
        burst = pd.read_csv(
            "_".join(basefile.split("_")[0:-1]) + "_burst.txt",
            infer_datetime_format=True,
        )

        burst = burst.rename(columns={"Time": "time"}).set_index("time")

        burst = burst.to_xarray()

        r = np.shape(burst.Pressure)[0]
        if r % ds.attrs["samples_per_burst"]:
            print(
                "Number of rows is not a multiple of samples_per_burst; truncating to last full burst"
            )
            mod = r % ds.attrs["samples_per_burst"]
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

        if "Wave" in dsburst:
            dsburst = dsburst.drop("Wave")

        ds = ds.sel(time=dsburst.sel(sample=0).time)
        # sample gets added with the .sel above
        ds = ds.drop("sample")
        # drop the burst average value to replae with burst data
        ds = ds.drop("P_1")
        ds = xr.merge([ds, dsburst])
        # dsburst.attrs = ds.attrs

    # Set burst interval, [sec], USER DEFINED in instrument attr
    # for continuous mode sampling
    if (ds.attrs["sample_mode"] == "CONTINUOUS") and ("wave_interval" in ds.attrs):
        ds.attrs["samples_per_burst"] = (
            ds.attrs["wave_interval"] / ds.attrs["sample_interval"]
        )
        ds.attrs["burst_interval"] = ds.attrs["wave_interval"]
        ds.attrs["burst_length"] = (
            ds.attrs["samples_per_burst"] * ds.attrs["sample_interval"]
        )
        r = np.shape(ds.P_1)[0]
        if r % ds.attrs["wave_interval"]:
            print(
                "Number of rows is not a multiple of wave_interval; truncating to last full burst"
            )
            mod = r % ds.attrs["wave_interval"]
            ds = ds.sel(time=ds.time[0:-mod])
        ds["timenew"] = xr.DataArray(
            ds.time[0 :: int(ds.attrs["wave_interval"])].values, dims="timenew"
        )
        ds["P_1"] = xr.DataArray(
            np.reshape(ds.P_1.values, (-1, int(ds.attrs["wave_interval"]))),
            dims=["timenew", "sample"],
        )
        ds = ds.rename({"time": "timeold"})
        ds = ds.rename({"timenew": "time"})
        ds = ds.drop("timeold")

    ds = utils.shift_time(ds, 0)

    # configure file
    cdf_filename = ds.attrs["filename"] + "-raw.cdf"

    ds.to_netcdf(cdf_filename, unlimited_dims=["time"])

    print("Finished writing data to %s" % cdf_filename)

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

    # print(ds[var])
    return ds
