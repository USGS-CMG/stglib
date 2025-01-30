import glob
import warnings

import numpy as np
import pandas as pd
import xarray as xr
import yaml
from tqdm import tqdm

from .aqd import aqdutils
from .core import filter, qaqc, utils, waves


def dat_to_cdf(metadata):
    """Read in Geolux Wave Radar data files and convert raw data into netcdf file"""

    basefile = metadata["basefile"]

    datfiles = glob.glob(f"{basefile}*.dat") + glob.glob(f"{basefile}*.txt")

    dflst = []

    print(f"reading in Geolux wave radar raw data files from {basefile}")
    for k in tqdm(datfiles):

        # Read in wave radar .dat or .txt files and make list of pandas dataframes.
        df = pd.read_csv(k, header=1, engine="pyarrow")
        df = df[2::]

        dflst.append(df)

    # Concatenate dataframes and create time index.
    df = pd.concat(dflst)
    df = df.reset_index(drop=True)
    df["time"] = pd.to_datetime(df["TIMESTAMP"], format="ISO8601")
    df.set_index("time", inplace=True)
    df.drop(columns="TIMESTAMP", inplace=True)

    for j in df.columns:
        if "RECORD" in j:
            df[j] = df[j].astype(int)
        else:
            df[j] = df[j].astype(float)

    # Convert concatenated dataframe to xarray data set and sort by time.
    ds = df.to_xarray()
    ds = ds.sortby("time")

    # write out metadata first, then deal exclusively with xarray attrs
    ds = utils.write_metadata(ds, metadata)

    del metadata

    ds = utils.ensure_cf(ds)

    ds["time"] = pd.DatetimeIndex(ds["time"])

    ds = utils.shift_time(ds, 0)
    # configure file
    cdf_filename = ds.attrs["filename"] + "-raw.cdf"

    ds.to_netcdf(cdf_filename, unlimited_dims=["time"])

    print("Finished writing data to %s" % cdf_filename)

    return ds


def cdf_to_nc(cdf_filename):
    """
    Load a raw .cdf file and generate a processed .nc file
    """

    # Load raw .cdf data
    ds = xr.open_dataset(cdf_filename)

    # Clip data to in/out water times or via good_ens
    ds = utils.clip_ds(ds)

    # make sure time index is unique
    idx = pd.Index(ds["time"])
    if len(idx) < len(ds["time"]):
        ds = ds.drop_duplicates(dim="time")

    # rename variables
    ds = ds_rename_vars(ds)

    # drop some vars after renaming
    for k in [
        "record",
        "avg_water_level",
    ]:
        if k in ds:
            ds = ds.drop_vars(k)

    ds = utils.create_z(ds)

    if "sample_rate" in ds.attrs and "sample_interval" not in ds.attrs:
        ds.attrs["sample_interval"] = 1 / ds.attrs["sample_rate"]

    elif "sample_interval" in ds.attrs and "sample_rate" not in ds.attrs:
        ds["sample_rate"] = 1 / ds.attrs["sample_interval"]

    else:  # find sample_rate and sample interval from data
        print("need to add find sample_rate and sample interval from time dim")
        ds.attrs["sample_interval"] = ds["time"].isel(time=slice(0, 10)).diff(
            dim="time"
        ).median().values / np.timedelta64(1, "s")
        ds.attrs["sample_rate"] = 1 / ds.attrs["sample_interval"]

    # fill bad values and gaps
    if "fill_value" in ds.attrs:  # if fill value is in data convert to Nans
        for var in ds.data_vars:
            ds[var] = ds[var].where(~(ds[var] == float(ds.attrs["fill_value"])), np.nan)

    for v in ["water_level", "brange"]:
        if v in ds:
            ds = qaqc.trim_min(ds, v)
            ds = qaqc.trim_max(ds, v)
            ds = qaqc.trim_min_diff(ds, v)
            ds = qaqc.trim_min_diff_pct(ds, v)
            ds = qaqc.trim_max_diff(ds, v)
            ds = qaqc.trim_max_diff_pct(ds, v)
            ds = qaqc.trim_bad_ens(ds, v)

    # fill gaps in time series
    ds = fill_time_gaps(ds)

    if "filtered_wl" in ds.attrs and ds.attrs["filtered_wl"].lower() == "true":
        ds = create_filtered_water_level(ds)

    ds = ds_add_attrs(ds)
    ds = utils.add_min_max(ds)
    ds = utils.add_start_stop_time(ds)
    # ds = utils.add_delta_t(ds)
    ds = utils.ds_add_lat_lon(ds)

    # Write to .nc file
    print("Writing cleaned/trimmed burst data and averaged burst data to .nc file")
    nc_filename = ds.attrs["filename"] + "cont-cal.nc"

    ds["time"] = ds["time"].dt.round("ms")

    ds.to_netcdf(
        nc_filename, unlimited_dims=["time"], encoding={"time": {"dtype": "double"}}
    )

    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])

    print("Done writing netCDF file", nc_filename)

    return ds


def nc_to_waves(nc_filename):
    """Calculate Wave Statistics for Geolux wave radar data"""
    ds = xr.open_dataset(nc_filename)

    # check to see if need to make wave burst from continuous data
    """if (ds.attrs["sample_mode"].lower() == "CONTINUOUS") and ("wave_interval" in ds.attrs):
        # check for wave_start_time attrs"""
    if "wave_start_time" in ds.attrs:
        print(
            f"trimming continuous data to start at users specified wave start time {ds.attrs['wave_start_time']}"
        )
        ds = ds.sel(
            time=slice(np.datetime64(ds.attrs["wave_start_time"]), ds["time"][-1])
        )

    # make elevation variable
    elev = ds.attrs["initial_instrument_height"] - ds["brange"]

    # set tolerance for filling gaps in wave burst data
    if "wavedat_tolerance" not in ds.attrs:
        ds.attrs["wavedat_tolerance"] = "20 s"

    tol = ds.attrs["wavedat_tolerance"]

    # replace NaNs within tolerance so can do wave stats
    elev = elev.where(~elev.isnull(), drop=True)
    elev = elev.reindex_like(ds["brange"], method="nearest", tolerance=tol)
    ds["elev"] = xr.DataArray(elev)

    # make wave burst ncfile from continuous data if wave_interval is specified
    ds = make_wave_bursts(ds)

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

    spec = waves.make_waves_ds_elev(ds)

    ds = ds_add_waves_history(ds)

    for k in ["wp_peak", "wh_4061", "wp_4060", "sspec"]:
        ds[k] = spec[k]

    # ds = utils.create_water_depth_var(ds)
    ds["water_depth"] = xr.DataArray(ds["elev"].squeeze().mean(dim="sample"))
    ds["water_depth"].attrs.update(
        {
            "long_name": "Total water depth",
            "units": "m",
            "standard_name": "sea_floor_depth_below_sea_surface",
            "epic_code": 3,
        }
    )

    ds["water_level"] = ds["water_level"].mean(dim="sample")
    ds["water_level"].attrs.update(
        {
            "units": "m",
            "long_name": "Water level NAVD88",
            "standard_name": "sea_surface_height_above_geopotential_datum",
            "geopotential_datum_name": "NAVD88",
        }
    )

    for k in [
        "water_level_filt",
        "sample",
        "brange",
        "elev",
    ]:
        if k in ds:
            ds = ds.drop_vars(k)

    # ds = qaqc.drop_vars(ds)

    # ds = drop_unused_dims(ds)

    ds = utils.trim_max_wp(ds)

    ds = utils.trim_min_wh(ds)

    ds = utils.trim_max_wh(ds)

    ds = utils.trim_wp_ratio(ds)

    # Add attrs
    ds = utils.ds_add_wave_attrs(ds)

    # clean-up attrs
    # ds = drop_attrs(ds)

    # assign min/max (need to do this after trimming):
    ds = utils.add_min_max(ds)

    # round time to minutes to make sure fits in dtype i4. will be fine for wave burst start times
    ds["time"] = ds["time"].dt.round("min")

    nc_filename = ds.attrs["filename"] + "s-a.nc"

    ds.to_netcdf(nc_filename, unlimited_dims=["time"])
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])

    print("Done writing netCDF file", nc_filename)

    return ds


def make_wave_bursts(ds):
    """Shape continuously sampled data into wave bursts using user specified interval"""
    # wave_interval is [sec] interval for wave statistics for continuous data
    ds.attrs["wave_samples_per_burst"] = int(
        ds.attrs["wave_interval"] / ds.attrs["sample_interval"]
    )

    r = np.shape(ds.brange)[0]
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

    for v in ["brange", "water_level", "water_level_filt", "elev"]:
        if v in ds:
            attrsbak = ds[v].attrs
            ds[v] = xr.DataArray(
                np.reshape(ds[v].values, (-1, int(ds.attrs["wave_samples_per_burst"]))),
                dims=["timenew", "samplenew"],
            )
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


def ds_rename_vars(ds):
    """Rename variables for CF output"""
    # modified from exo.ds_rename_vars
    varnames = {
        "RECORD": "record",
        "Level_10Hz": "water_level",
        "Avg_Level_10Hz": "avg_water_level",
        "Distance_10Hz": "brange",
    }

    # check to make sure they exist before trying to rename
    newvars = {}
    for k in varnames:
        if k in ds:
            newvars[k] = varnames[k]

    return ds.rename(newvars)


def ds_add_attrs(ds):
    """Add variable attribute to some Geolux specific variables"""
    # modified from exo.ds_add_attrs
    ds = utils.ds_coord_no_fillvalue(ds)

    ds["time"].attrs.update(
        {"standard_name": "time", "axis": "T", "long_name": "time (UTC)"}
    )

    # ds["sample"].attrs.update({"units": "1", "long_name": "Sample in burst"})

    ds["water_level"].attrs.update(
        {
            "units": "m",
            "long_name": "Water level NAVD88",
            "standard_name": "sea_surface_height_above_geopotential_datum",
            "geopotential_datum_name": "NAVD88",
        }
    )

    ds["brange"].attrs.update(
        {
            "units": "m",
            "long_name": "sensor range to boundary",
            "standard_name": "altimeter_range",
        }
    )

    return ds


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


def create_filtered_water_level(ds):
    """Create filetered water level, fill missing water level values temporarily to generate filter water level if largest gaps are less than tolerance"""

    print("Add filtered water level variable")
    # Nans in water level variable, to calculate filetered water level need to try to fill Nans if gaps aren't too big
    nnans = ds["water_level"].isnull().sum().values

    # set tolerance for filling gaps in water level data for calculating filtered water level
    if "wlfilt_tolerance" not in ds.attrs:
        ds.attrs["wlfilt_tolerance"] = "60 s"

    tol = ds.attrs["wlfilt_tolerance"]

    # for water_level temporarily replace Nans with interpolated vals to find filtered water level then replace back
    dsf = ds.where(~ds["water_level"].isnull(), drop=True)
    dsf = dsf.reindex_like(ds, method="nearest", tolerance=tol)

    dsf = utils.create_filtered_water_level_var(dsf)

    # fill with Nans where water_level is Nan.
    ds["water_level_filt"] = dsf["water_level_filt"]
    ds["water_level_filt"] = ds["water_level_filt"].where(
        ~ds["water_level"].isnull(), np.nan
    )

    if dsf["water_level"].isnull().sum() > 0:
        fwl_note = "Gaps in valid water level exceed tolerance and filtered water level will be null values"
    else:
        fwl_note = f"To calculate filtered water level {nnans} out of {len(ds['water_level'].values)} water level values were temporarily filled with nearest neighbor using tolerance = {tol}s"

    print(fwl_note)

    if "note" in ds["water_level_filt"].attrs:
        notetxt = ds["water_level_filt"].attrs["note"] + "; " + fwl_note
    else:
        notetxt = fwl_note

    ds["water_level_filt"].attrs.update({"note": notetxt})

    return ds


def ds_add_waves_history(ds):
    histtext = (
        "Wave statistics computed using scipy.signal.welch(), "
        f"Since this instrument ({ds.attrs['instrument_type']}) measures the sea-surface elevation directly no cutoff frequecncy is assigned, but a maximum frequency of the lesser of Nyquist (sample_rate/2) or 2 Hz is used."
    )

    return utils.insert_history(ds, histtext)
