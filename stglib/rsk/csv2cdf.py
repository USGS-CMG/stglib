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

    ds.attrs["serial_number"] = meta["instrument"]['serial']

    for k in ds.data_vars:
        if " " in k:
            ds = ds.rename({k: k.replace(" ", "_")})

    if "Turbidity" in ds:
        ds = ds.rename({"Turbidity": "Turb"})

        ds["Turb"].attrs.update(
            {"units": "Nephelometric turbidity units (NTU)", "long_name": "Turbidity"}
        )

    if "Pressure"in ds:
        ds = ds.rename({"Pressure": "P_1"})

    if "Sea_pressure" in ds:
        ds = ds.drop("Sea_pressure")

    if "Depth" in ds:
        ds = ds.drop("Depth")

    # Set burst interval, [sec], USER DEFINED in instrument attr
    # for continuous mode sampling
    if "wave_interval" in ds.attrs:
        print("in wave interval code")
        ds.attrs["sample_interval"] = int(meta["sampling"]["period"]) / 1000
        ds.attrs["sample_mode"] = meta["sampling"]["mode"]
        ds.attrs["samples_per_burst"] = ds.attrs["wave_interval"] / ds.attrs["sample_interval"]
        ds.attrs["burst_interval"] = ds.attrs["wave_interval"]
        ds.attrs["burst_length"] = (
            ds.attrs["samples_per_burst"] * ds.attrs["sample_interval"]
        )
        r = np.shape(ds.P_1)[0]
        if r % ds.attrs["wave_interval"]:
            print("Number of rows is not a multiple of wave_interval; truncating to last full burst")
            mod = r % ds.attrs["wave_interval"]
            ds = ds.sel(
                time=ds.time[0 : -mod]
            )
        ds["timenew"] = xr.DataArray(ds.time[0::int(ds.attrs["wave_interval"])].values, dims='timenew')
        ds["P_1"] = xr.DataArray(np.reshape(ds.P_1.values, (-1, int(ds.attrs["wave_interval"]))), dims=['timenew', 'sample'])
        ds = ds.rename({"time": "timeold"})
        ds = ds.rename({"timenew": "time"})
        ds = ds.drop("timeold")

    ds = utils.shift_time(ds, 0)

    # configure file
    cdf_filename = ds.attrs["filename"] + "-raw.cdf"

    ds.to_netcdf(cdf_filename, unlimited_dims=["time"])

    print("Finished writing data to %s" % cdf_filename)

    return ds
