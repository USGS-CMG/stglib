import math

import numpy as np
import pandas as pd
import xarray as xr

from ..core import utils
from . import sgutils


def cdf_to_nc(cdf_filename, atmpres=None):
    """
    Load a raw .cdf file and generate a processed .nc file
    """

    # Load raw .cdf data
    ds = xr.open_dataset(cdf_filename)

    # Atmospheric pressure correction
    if atmpres is not None:
        ds = sgutils.atmos_correct_burst(ds, atmpres)

    # If .wb file, split into tide bursts
    if ds.attrs["file_type"] == ".wb":
        ds = avg_tide_bursts(ds)

    # Rename variables to CF compliant names
    ds = ds_rename_vars(ds)

    # Add attributes
    ds = sgutils.ds_add_attrs(ds)

    # Edit metadata depending on file type
    if ds.attrs["file_type"] == ".tid":
        ds = ds_drop_tid(ds)
    elif ds.attrs["file_type"] == ".wb":
        ds = ds_drop_wb(ds)

    # Call QAQC
    ds = sgutils.sg_qaqc(ds)

    # Run utilities
    ds = utils.clip_ds(ds)
    ds = utils.ds_add_lat_lon(ds)
    ds = utils.create_nominal_instrument_depth(ds)
    ds = utils.create_z(ds)
    ds = utils.add_start_stop_time(ds)
    ds = utils.add_min_max(ds)
    ds = utils.add_delta_t(ds)

    # Write to .nc file
    print("Writing cleaned/trimmed data to .nc file")
    nc_filename = ds.attrs["filename"] + "s-tide-a.nc"

    ds.to_netcdf(
        nc_filename, unlimited_dims=["time"], encoding={"time": {"dtype": "i4"}}
    )
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])

    print(f"Done writing netCDF file {nc_filename}")


def ds_rename_vars(ds):
    """
    Rename variables to be CF compliant
    """
    varnames = {"Temp": "T_28"}

    # Check to make sure they exist before trying to rename
    newvars = {}
    for k in varnames:
        if k in ds:
            newvars[k] = varnames[k]
    return ds.rename(newvars)


def ds_drop_tid(ds):
    """
    Drop global attribute metadata not needed for .tid file
    """
    gatts = [
        "WaveInterval",
        "WaveIntervalUnits",
        "WaveSamples",
        "sample_rate",
        "sample_rate_units",
        "BurstDuration",
        "BurstDurationUnits",
        "WaveBurstsPerDay",
        "NumberOfWaveBursts",
        # "calculated_wave_interval",
        # "calculated_wave_interval_units",
    ]

    # Check to make sure they exist
    for k in gatts:
        if k in ds.attrs:
            del ds.attrs[k]
    return ds


def ds_drop_wb(ds):
    """
    Drop global attribute metadata not needed for .wb file
    """
    gatts = [
        "WaveInterval",
        "WaveIntervalUnits",
        "WaveSamples",
        "sample_rate",
        "sample_rate_units",
        "BurstDuration",
        "BurstDurationUnits",
        "WaveBurstsPerDay",
        "NumberOfWaveBursts",
        # "calculated_wave_interval",
        # "calculated_wave_interval_units",
        "TideInterval",
        "TideIntervalUnits",
        "TideDuration",
        "TideDurationUnits",
        "TideSamplesPerDay",
        "NumberOfTideMeasurements",
    ]

    # Check to make sure they exist
    for k in gatts:
        if k in ds.attrs:
            del ds.attrs[k]
    return ds


def avg_tide_bursts(ds):

    # Calculate how many rows to subdivide each wave burst (round down to trim burst)
    rows = math.floor(
        float(ds.attrs["BurstDuration"]) / float(ds.attrs["calculated_tide_interval"])
    )

    # Calculate how many columns to subdivide each wave burst
    cols = int(
        float(ds.attrs["calculated_tide_interval"]) * float(ds.attrs["sample_rate"])
    )

    # Calculate number of values that can be used from each wave burst
    no_values = int(
        float(ds.attrs["calculated_tide_interval"])
        * float(ds.attrs["sample_rate"])
        * rows
    )

    # Calculate number of values to average based on tide duration
    values_avg = int(
        float(ds.attrs["calculated_tide_duration"]) * float(ds.attrs["sample_rate"])
    )

    # Define time interval
    delta_t = int(ds.attrs["calculated_tide_interval"])
    delta_t = str(delta_t) + "s"

    # Reshape P_1 and P_1ac pressure bursts and average
    P1_avg = []
    P1ac_avg = []
    timestamp = []
    for t in range(0, ds["time"].size):

        # Make timestamp
        new_time = np.array(
            pd.date_range(start=ds.time[t].values, periods=rows, freq=delta_t)
        )
        timestamp.append(new_time)

        # Trim incomplete bursts
        P1_trim_burst = ds["P_1"].isel(time=t, sample=slice(0, no_values))
        P1ac_trim_burst = ds["P_1ac"].isel(time=t, sample=slice(0, no_values))

        # Reshape
        new_P1_burst = np.reshape(P1_trim_burst.values, (rows, cols))
        new_P1ac_burst = np.reshape(P1ac_trim_burst.values, (rows, cols))

        # Average burst for tide duration
        for j in range(0, new_P1_burst.shape[0]):
            avg_pres_P1 = np.mean(new_P1_burst[j, slice(0, values_avg)])
            avg_pres_P1ac = np.mean(new_P1ac_burst[j, slice(0, values_avg)])
            P1_avg.append(avg_pres_P1)
            P1ac_avg.append(avg_pres_P1ac)

    # Append all timestamps to array
    date = np.concatenate(timestamp)

    # Copy global attributes
    attr_copy = ds.attrs.copy()

    # Delete old xarray because no longer need burst data
    del ds

    # Make new xarray with averaged data
    ds = xr.Dataset(
        data_vars=dict(P_1=(["time"], P1_avg), P_1ac=(["time"], P1ac_avg)),
        coords=dict(time=date),
    )

    # Add back global attributes
    ds.attrs = attr_copy

    return ds
