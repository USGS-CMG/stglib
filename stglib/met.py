import numpy as np
import pandas as pd
import xarray as xr
from scipy import stats

from .core import attrs, qaqc, utils


def read_mayfly(filnam, skiprows=7):
    """Read csv file from EnviroDIY Mayfly datalogger into an xarray
    Dataset.

    skiprows : int, optional
        How many header rows to skip. Default 7.
    """

    df = pd.read_csv(
        filnam,
        skiprows=skiprows,
        header=0,
        na_values=[-9999],
    )
    df["Date and Time in UTC"] = pd.to_datetime(df["Date and Time in UTC"])
    df.rename(columns={"Date and Time in UTC": "time"}, inplace=True)
    df.set_index("time", inplace=True)

    return df.to_xarray()


def read_campbell(filnam):
    """Read dat file from Campbell Scientific datalogger into an xarray
    Dataset.
    """

    df = pd.read_csv(
        filnam,
        header=1,
        skiprows=[2, 3],
        na_values=["NAN"],
    )
    df["TIMESTAMP"] = pd.to_datetime(df["TIMESTAMP"])
    df.rename(columns={"TIMESTAMP": "time"}, inplace=True)
    df.set_index("time", inplace=True)

    return df.to_xarray()


# Make raw CDF
def csv_to_cdf(metadata):
    basefile = metadata["basefile"]

    if metadata["datalogger"].lower() == "mayfly":
        ds = read_mayfly(basefile + ".csv", skiprows=metadata["skiprows"])
        metadata.pop("skiprows")
    elif metadata["datalogger"].lower() == "campbell":
        ds = read_campbell(basefile + ".dat")

    ds = utils.write_metadata(ds, metadata)

    ds = utils.ensure_cf(ds)

    # configure file
    cdf_filename = ds.attrs["filename"] + "-raw.cdf"

    ds.to_netcdf(cdf_filename, unlimited_dims=["time"])

    print("Finished writing data to %s" % cdf_filename)

    return ds


# Process data and write to .nc file
def cdf_to_nc(cdf_filename):
    """
    Load a raw .cdf file and generate a processed .nc file
    """

    ds = xr.open_dataset(cdf_filename)

    # remove units in case we change and we can use larger time steps
    ds.time.encoding.pop("units")

    # Drop unneeded variables
    ds = met_drop_vars(ds)

    # Rename variables to CF compliant names
    ds = ds_rename_vars(ds)

    # Convert vars from int64 to float (doesn't touch coords like time)
    ds = ds.astype(float)

    # Fill time gaps
    ds.attrs["sample_interval"] = stats.mode(
        np.diff(ds.time) / np.timedelta64(1, "s")
    ).mode

    ds.attrs["sample_rate"] = 1 / ds.attrs["sample_interval"]  # Hz

    ds = fill_time_gaps(ds)

    # Remove bad rows. Needs to happen before direction corrections.
    ds = qaqc.call_qaqc(ds)

    # Apply direction offset and magnetic declination correction
    wind_vars = [
        "WD_min",
        "WD_410",
        "WD_gust",
        "wind_dir",
    ]

    for var in wind_vars:
        if var in ds:

            # If sensor wasn't pointing to magnetic north, apply offset to direction
            if "dir_offset" in ds.attrs:

                ds[var] = ds[var] + ds.attrs["dir_offset"].astype(float)
                ds = utils.insert_history(
                    ds, f"Applied dir_offset of {ds.attrs['dir_offset']}"
                )

            # Convert direction from magnetic to true with magnetic declination
            ds[var] = ds[var] + ds.attrs["magnetic_variation"].astype(float)
            ds[var] = ds[var].round(0)
            ds[var] = ds[var] % 360

            ds = utils.insert_history(
                ds,
                f"Rotated directions from magnetic north to true north by applying magnetic_variation of {ds.attrs['magnetic_variation']}",
            )

    # Run utilities
    ds = utils.create_z(ds)
    ds = utils.clip_ds(ds)
    ds = utils.add_start_stop_time(ds)
    ds = utils.ds_add_lat_lon(ds)
    ds = utils.add_min_max(ds)
    ds = utils.add_delta_t(ds)

    # Add attributes
    ds = attrs.ds_add_attrs(ds)
    ds = ds_add_var_attrs(ds)

    # Write to .nc file
    print("Writing cleaned/trimmed data to .nc file")
    nc_filename = ds.attrs["filename"] + "-a.nc"

    ds.to_netcdf(
        nc_filename, unlimited_dims=["time"], encoding={"time": {"dtype": "i4"}}
    )
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])
    print("Done writing netCDF file", nc_filename)


def fill_time_gaps(ds):
    """Fill any gaps in time-series using to make time even"""

    sr = ds.attrs["sample_rate"]
    sims = 1 / sr * 1000
    pds = (
        int(
            (ds["time"][-1].values - ds["time"][0].values)
            / (sims * np.timedelta64(1, "ms"))
        )
        + 1
    )
    idx = pd.date_range(str(ds["time"][0].values), periods=pds, freq=f"{sims}ms")

    # make sure time index is unique
    ds = ds.drop_duplicates(dim="time")

    ds = ds.reindex(time=idx)

    return ds


def ds_rename_vars(ds):

    wxt_names = {
        # Mayfly + WXT
        "WXTDn": "WD_min",
        "WXTDm": "WD_410",
        "WXTDx": "WD_gust",
        "WXTSn": "WS_min",
        "WXTSm": "WS_401",
        "WXTSx": "WG_402",
        "WXTTa": "T_21",
        "WXTUa": "RH_910",
        "WXTPa": "BPR_915",
        "WXTRc": "Rn_963",
        # Campbell + WXT
        "WindDir_lull": "WD_min",
        "WindDir_avg": "WD_410",
        "WindDir_gust": "WD_gust",
        "WindSpeed_lull": "WS_min",
        "WindSpeed_avg": "WS_401",
        "WindSpeed_gust": "WG_402",
        "Temp": "T_21",
        "RH": "RH_910",
        "Baro": "BPR_915",
        "R_amt": "Rn_963",
        "R_dur": "rain_duration",
        "R_int": "rain_rate",
        "H_amt": "hail_amount",
        "H_dur": "hail_duration",
        "H_int": "hail_rate",
    }

    clivue_names = {
        # Campbell + ClimaVue
        "PTemp_C_Avg": "internal_temp",
        "SlrFD_W": "solar_flux_density",
        "Rain_mm_Tot": "rain_amount",
        "Strikes_Tot": "light_strikes",
        "Dist_km": "strike_distance",
        "WS_ms": "wind_speed",
        "WindDir": "wind_dir",
        "MaxWS_ms": "wind_gust",
        "AirT_C": "air_temp",
        "VP_mbar": "vapor_pressure",
        "BP_mbar": "baro_pressure",
        "RH": "relative_humidity",
        "RHT_C": "humidity_sensor_temp",
        "TiltNS_deg": "tilt_NS",
        "TiltWE_deg": "tilt_WE",
        "SlrTF_MJ_Tot": "solar_total_flux",
        "Invalid_Wind": "wind_error",
    }

    # Check to make sure vars exist before trying to rename
    newvars = {}

    if ds.attrs["instrument_type"].lower() == "wxt":
        for k in wxt_names:
            if k in ds:
                newvars[k] = wxt_names[k]
    elif ds.attrs["instrument_type"].lower() == "climavue":
        for k in clivue_names:
            if k in ds:
                newvars[k] = clivue_names[k]

    return ds.rename(newvars)


def met_drop_vars(ds):

    # Drop unneeded variables
    var_list = [
        "SampNum",
        "Battery",
        "BoardTemp",
        "signalPercent",
        "RECORD",
        "panel_temp",
        "power_in",
        "lithium_battery",
        "memory_free",
        "BattV_Max",
        "CVMeta",
    ]

    # Will ignore errors if variable is not in dataset
    ds = ds.drop_vars(var_list, errors="ignore")

    return ds


# Add initial height to specific variables
def ds_add_var_attrs(ds):

    for name in ds.variables:

        if (name not in ds.coords) and ("time" not in name):
            # don't include for coordinates that are also variables

            var = ds[name]

            var.attrs["initial_instrument_height"] = ds.attrs[
                "initial_instrument_height"
            ]
            var.attrs["height_depth_units"] = "m"
            if "initial_instrument_height_note" in ds.attrs:
                var.attrs["initial_instrument_height_note"] = ds.attrs[
                    "initial_instrument_height_note"
                ]

    return ds
