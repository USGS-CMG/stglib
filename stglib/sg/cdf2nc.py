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

    # If .wb file, split into tide bursts
    if ds.attrs["file_type"] == ".wb":
        ds = avg_tide_bursts(ds)

    # Atmospheric pressure correction
    if atmpres is not None:
        ds = utils.atmos_correct(ds, atmpres)

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


# def ds_drop_vars(ds):
#     """
#     Drop unneeded variables
#     """
#     varnames = ["sample", "burst_number","P_1"]

#     # Check to make sure they exist
#     for k in varnames:
#         if k in ds:
#             ds = ds.drop[k]
#     return ds


def ds_drop_tid(ds):
    """
    Drop global attribute metadata not needed for .tid file
    """
    gatts = [
        "WaveInterval",
        "WaveIntervalUnits",
        "WaveSamples",
        "WaveSampleRate",
        "WaveSampleRateUnits",
        "BurstDuration",
        "BurstDurationUnits",
        "WaveBurstsPerDay",
        "NumberOfWaveBursts",
        "calculated_wave_interval",
        "calculated_wave_interval_units",
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
        "WaveSampleRate",
        "WaveSampleRateUnits",
        "BurstDuration",
        "BurstDurationUnits",
        "WaveBurstsPerDay",
        "NumberOfWaveBursts",
        "calculated_wave_interval",
        "calculated_wave_interval_units",
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
        float(ds.attrs["calculated_tide_interval"]) * float(ds.attrs["WaveSampleRate"])
    )

    # Calculate number of values that can be used from each wave burst
    no_values = int(
        float(ds.attrs["calculated_tide_interval"])
        * float(ds.attrs["WaveSampleRate"])
        * rows
    )

    # Calculate number of values to average based on tide duration
    values_avg = int(
        float(ds.attrs["calculated_tide_duration"]) * float(ds.attrs["WaveSampleRate"])
    )

    # Define time interval
    delta_t = int(ds.attrs["calculated_tide_interval"])
    delta_t = str(delta_t) + "s"

    # Reshape pressure bursts and average
    burst = []
    timestamp = []
    for t in range(0, ds["time"].size):

        # Make timestamp
        new_time = np.array(
            pd.date_range(start=ds.time[t].values, periods=rows, freq=delta_t)
        )
        timestamp.append(new_time)

        # Trim incomplete bursts
        indiv_burst = ds["P_1"].isel(time=t, sample=slice(0, no_values))

        # Reshape
        new_burst = np.reshape(indiv_burst.values, (rows, cols))

        # Average burst for tide duration
        # for j in range(0,new_burst.shape[0]):
        #     avg_pres=np.mean(new_burst.isel(time=j, sample=slice(0,values_avg))))
        #     #avg_pres=np.mean(np.array(ds["P_1"].isel(time=j, sample=slice(0,values_avg))))
        #     # avg_pres = avg_pres.tolist()
        #     burst.append(avg_pres)

        # Average burst for tide duration
        for j in range(0, new_burst.shape[0]):
            avg_pres = np.mean(new_burst[j, slice(0, values_avg)])
            burst.append(avg_pres)

    # Append all timestamps to array
    date = np.concatenate(timestamp)

    # Copy global attributes
    attr_copy = ds.attrs.copy()

    # Delete old xarray because no longer need burst data
    del ds

    # Make new xarray with averaged data
    ds = xr.Dataset(
        data_vars=dict(P_1=(["time"], burst)),
        coords=dict(time=date),
    )

    # Add back global attributes
    ds.attrs = attr_copy

    return ds
