from __future__ import division, print_function

import warnings

import numpy as np
import pandas as pd
import xarray as xr

from .core import utils

# Read WXT
def read_wxt(filnam, skiprows=7, encoding="utf-8"):
    """Read data from a Vaisala WXT met .csv file into an xarray
    Dataset.
    Parameters
    ----------
    filnam : string
        The filename
    skiprows : int, optional
        How many header rows to skip. Default 7
    encoding : string, optional
        File encoding. Default 'utf-8'
    Returns
    -------
    xarray.Dataset
        An xarray Dataset of the WXT data
    """

    df = pd.read_csv(
        filnam,
        skiprows=skiprows,
        header=0,
        na_values=[-9999],
        encoding=encoding,
        index_col=False,
    )
    df["Date and Time in UTC"] = pd.to_datetime(df["Date and Time in UTC"])
    df.rename(columns={"Date and Time in UTC": "time"}, inplace=True)
    df.index.names = ["time"]
    wxt = xr.Dataset.from_dataframe(df)
    return wxt


# Make raw CDF
def csv_to_cdf(metadata):

    basefile = metadata["basefile"]

    ds = read_wxt(basefile + ".csv", skiprows=metadata["skiprows"])
    metadata.pop("skiprows")
    ds = utils.write_metadata(ds, metadata)
    # ds['time'] = xr.DataArray(time, dims='time')

    # configure file
    cdf_filename = ds.attrs["filename"] + "-raw.cdf"

    ds.to_netcdf(cdf_filename, unlimited_dims=["time"])

    print("Finished writing data to %s" % cdf_filename)

    return ds


# Process data and write to .nc file
def cdf_to_nc(cdf_filename):
    """
    Load a "raw" .cdf file and generate a processed .nc file
    """

    # Load raw .cdf data
    ds = xr.open_dataset(cdf_filename)

    # Get rid of unneeded variables
    for k in [
        "SampNum",
        "Battery",
        "BoardTemp",
        "signalPercent",
    ]:
        if k in ds:
            ds = ds.drop_vars(k)

    # Add height as variable
    ds["height"] = xr.DataArray(
        [ds.attrs["initial_instrument_height"]],
        dims=("height"),
        name="height",
        attrs={
            "units": "m",
            "long_name": "Measurement Elevation",
            "standard_name": "height",
            "positive": "up",
            "axis": "Z",
            "epic_code": "18",
        },
    )

    # Rename variables to CF compliant names
    ds = ds_rename_vars(ds)

    # Convert data types to float 32
    ds = ds_convertfloat(ds)

    # If sensor wasn't pointing to magnetic north, apply offset to direction
    if "dir_offset" in ds.attrs:
        if "WD_min" in ds:
            ds["WD_min"].values = ds["WD_min"].values + ds.attrs["dir_offset"]
        if "WD_410" in ds:
            ds["WD_410"].values = ds["WD_410"].values + ds.attrs["dir_offset"]
        if "WD_gust" in ds:
            ds["WD_gust"].values = ds["WD_gust"].values + ds.attrs["dir_offset"]

    # Convert direction from magnetic to true with magenetic declination
    if "WD_min" in ds:
        ds["WD_min"].values = ds["WD_min"].values + ds.attrs["magnetic_variation"]
        ds["WD_min"].values[ds["WD_min"].values < 0.0] = (
            ds["WD_min"].values[ds["WD_min"].values < 0.0] + 360.0
        )
        ds["WD_min"].values[ds["WD_min"].values > 360.0] = (
            ds["WD_min"].values[ds["WD_min"].values > 360.0] - 360.0
        )

    if "WD_410" in ds:
        ds["WD_410"].values = ds["WD_410"].values + ds.attrs["magnetic_variation"]
        ds["WD_410"].values[ds["WD_410"].values < 0.0] = (
            ds["WD_410"].values[ds["WD_410"].values < 0.0] + 360.0
        )
        ds["WD_410"].values[ds["WD_410"].values > 360.0] = (
            ds["WD_410"].values[ds["WD_410"].values > 360.0] - 360.0
        )

    if "WD_gust" in ds:
        ds["WD_gust"].values = ds["WD_gust"].values + ds.attrs["magnetic_variation"]
        ds["WD_gust"].values[ds["WD_gust"].values < 0.0] = (
            ds["WD_gust"].values[ds["WD_gust"].values < 0.0] + 360.0
        )
        ds["WD_gust"].values[ds["WD_gust"].values > 360.0] = (
            ds["WD_gust"].values[ds["WD_gust"].values > 360.0] - 360.0
        )

    # Run utilities
    ds = utils.add_start_stop_time(ds)
    ds = utils.ds_add_lat_lon(ds)
    ds = utils.add_min_max(ds)
    ds = utils.add_delta_t(ds)

    # Add attributes
    ds = ds_add_attrs(ds)

    # Write to .nc file
    print("Writing cleaned/trimmed data to .nc file")
    nc_filename = ds.attrs["filename"] + "-a.nc"

    ds.to_netcdf(
        nc_filename, unlimited_dims=["time"], encoding={"time": {"dtype": "i4"}}
    )
    utils.check_compliance(nc_filename)
    print("Done writing netCDF file", nc_filename)


# Rename variables to be CF compliant
def ds_rename_vars(ds):
    varnames = {
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
    }

    # Check to make sure they exist before trying to rename
    newvars = {}
    for k in varnames:
        if k in ds:
            newvars[k] = varnames[k]
    return ds.rename(newvars)


# Convert data from float 64 to float32
def ds_convertfloat(ds):
    for var in ds.variables:
        if ds[var].name != "time":
            ds[var] = ds[var].astype("float32")
    return ds


# Add attributes: units, standard name from CF website, epic code
def ds_add_attrs(ds):

    ds = utils.ds_coord_no_fillvalue(ds)

    ds["time"].attrs.update({"standard_name": "time", "axis": "T"})

    if "WD_min" in ds:
        ds["WD_min"].attrs.update(
            {"units": "degrees", "long_name": "minimum wind from direction"}
        )

    if "WD_410" in ds:
        ds["WD_410"].attrs.update(
            {
                "units": "degrees",
                "long_name": "mean wind from direction",
                "standard_name": "wind_from_direction",
                "epic_code": "410",
            }
        )

    if "WD_gust" in ds:
        ds["WD_gust"].attrs.update(
            {
                "units": "degrees",
                "long_name": "maximum wind from direction",
                "standard_name": "wind_gust_from_direction",
            }
        )

    if "WS_min" in ds:
        ds["WS_min"].attrs.update({"units": "m/s", "long_name": "minimum wind speed"})

    if "WS_401" in ds:
        ds["WS_401"].attrs.update(
            {
                "units": "m/s",
                "long_name": "mean wind speed",
                "standard_name": "wind_speed",
                "epic_code": "401",
            }
        )

    if "WG_402" in ds:
        ds["WG_402"].attrs.update(
            {
                "units": "m/s",
                "long_name": "maximum wind speed",
                "standard_name": "wind_speed_of_gust",
                "epic_code": "402",
            }
        )

    if "T_21" in ds:
        ds["T_21"].attrs.update(
            {"units": "degree_C", "standard_name": "air_temperature", "epic_code": "21"}
        )

    if "RH_910" in ds:
        ds["RH_910"].attrs.update(
            {"units": "%", "standard_name": "relative_humidity", "epic_code": "910"}
        )

    if "BPR_915" in ds:
        ds["BPR_915"].attrs.update(
            {"units": "pascals", "standard_name": "air_pressure", "epic_code": "915"}
        )

    if "Rn_963" in ds:
        ds["Rn_963"].attrs.update(
            {
                "units": "mm",
                "standard_name": "thickness_of_rainfall_amount",
                "epic_code": "963",
            }
        )

    #     add initial height information and fill values to variables
    def add_attributes(var, dsattrs):
        var.attrs.update(
            {
                "initial_instrument_height": dsattrs["initial_instrument_height"],
                "height_depth_units": "m",
                "initial_instrument_height_note": dsattrs[
                    "initial_instrument_height_note"
                ],
                "sensor_type": "Vaisala WXT536",
            }
        )

    for var in ds.variables:
        if ds[var].dtype == "float32":
            ds[var].encoding["_FillValue"] = 1e35
        elif ds[var].dtype == "int32":
            ds[var].encoding["_FillValue"] = -2147483648

    # don't include all attributes for coordinates that are also variables
    for var in ds.variables:
        if (var not in ds.coords) and ("time" not in var):
            add_attributes(ds[var], ds.attrs)

    return ds
