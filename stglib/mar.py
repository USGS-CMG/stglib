import math

import numpy as np
import pandas as pd
import xarray as xr
from scipy import stats

from .aqd import aqdutils
from .core import qaqc, utils


def read_mar(filnam, encoding="utf-8"):
    """Read data from a Marotte Tilt Current Meter csv file into an xarray
    Dataset.
    """
    df = pd.read_csv(
        filnam,
        skiprows=1,
        header=None,
        names=[
            "time",
            "speed",
            "heading",
            "speed_upper",
            "speed_lower",
            "tilt",
            "direction",
            "batt",
            "temp",
        ],
        encoding=encoding,
        index_col=False,
    )

    df["time"] = pd.to_datetime(df["time"])
    df = df.set_index("time")

    return df.to_xarray()


def csv_to_cdf(metadata):
    """
    Load a raw .csv file and generate a .cdf file
    """
    basefile = metadata["basefile"]

    ds = read_mar(basefile + ".csv")

    ds = utils.write_metadata(ds, metadata)

    ds = utils.ensure_cf(ds)

    # Configure file
    cdf_filename = ds.attrs["filename"] + "-raw.cdf"

    ds.to_netcdf(cdf_filename, unlimited_dims=["time"])

    print(f"Finished writing data to {cdf_filename}")

    return ds


def cdf_to_nc(cdf_filename):
    """
    Load a raw .cdf file and generate a processed .nc file
    """

    # Load raw .cdf data
    ds = xr.open_dataset(cdf_filename)

    # Remove units in case we change and we can use larger time steps
    ds.time.encoding.pop("units")

    # Calculate sample interval/rate from most frequent time difference between samples in case there are any time skips
    ds.attrs["sample_interval"] = stats.mode(np.diff(ds.time.values))[
        0
    ] / np.timedelta64(1, "s")
    ds.attrs["sample_rate"] = 1 / ds.attrs["sample_interval"]  # Hz

    # Calculate u and v velocities
    ds["u"], ds["v"] = utils.spd2uv(ds["speed"], ds["heading"])

    # Drop variables
    ds = ds.drop_vars(["speed_upper", "speed_lower", "direction", "tilt", "batt"])

    # Rename variables to CF compliant names
    ds = ds_rename_vars(ds)

    # Apply magnetic variation after renaming variables
    if ds.attrs["correct_mag_var"].upper() == "TRUE":
        ds = aqdutils.magvar_correct(ds)

    # If burst data, reshape into bursts
    if ds.attrs["sample_mode"].upper() == "BURST":
        ds = reshape_burst(ds)

    # Add attributes after magnetic variation correction
    ds = ds_add_attrs(ds)

    # QAQC
    ds = qaqc.call_qaqc(ds)

    # Run utilities
    ds = utils.create_z(ds)
    ds = utils.clip_ds(ds)
    ds = utils.ds_add_lat_lon(ds)
    ds = utils.create_nominal_instrument_depth(ds)
    ds = utils.add_start_stop_time(ds)
    ds = utils.add_min_max(ds)
    ds = utils.ds_coord_no_fillvalue(ds)

    # Write to .nc file
    print("Writing cleaned/trimmed data to .nc file")
    nc_filename = ds.attrs["filename"] + "b.nc"

    ds["time"].encoding["dtype"] = "i4"
    ds.to_netcdf(nc_filename, unlimited_dims=["time"])
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])

    print(f"Done writing netCDF file {nc_filename}")

    # Average continuous or burst file
    ds_avg = None

    if (
        ds.attrs["sample_mode"].upper() == "CONTINUOUS"
    ) and "average_interval" in ds.attrs:

        # Check for time gaps and fill before averaging
        ds_fill = fill_time_gaps(ds)

        # Average continuous data
        ds_avg = avg_cont(ds_fill)

    elif (
        ds.attrs["sample_mode"].upper() == "BURST"
        and ds.attrs["average_duration"] not in ds.attrs
    ):

        # Average entire burst
        ds_avg = ds.mean(dim="sample")

    elif (
        ds.attrs["sample_mode"].upper() == "BURST"
        and ds.attrs["average_duration"] in ds.attrs
    ):

        # Get number of samples to average per burst
        samp_per_burst = int(ds.attrs["average_duration"] * ds.attrs["sample_rate"])

        # Grab only number of samples needed
        ds_slice = ds.isel(sample=slice(0, samp_per_burst))

        # Average duration of burst
        ds_avg = ds_slice.mean(dim="sample")

    if ds_avg is not None:

        # Add back global attributes
        ds_avg.attrs = ds.attrs

        # Re-calculate average speed and heading from averaged u and v
        ds_avg["CS_300"], ds_avg["CD_310"] = utils.uv2spd(
            ds_avg["u_1205"], ds_avg["v_1206"]
        )

        # Re-run utils
        ds_avg = ds_add_attrs(ds_avg)
        ds_avg = utils.create_z(ds_avg)
        ds_avg = utils.ds_add_lat_lon(ds_avg)
        ds_avg = utils.create_nominal_instrument_depth(ds_avg)
        ds_avg = utils.add_start_stop_time(ds_avg)
        ds_avg = utils.add_min_max(ds_avg)
        ds_avg = utils.add_delta_t(ds_avg)
        ds_avg = utils.ds_coord_no_fillvalue(ds_avg)

        # Write to .nc file
        print("Writing averaged data to .nc file")
        nc_filename = ds_avg.attrs["filename"] + "b-a.nc"

        ds_avg["time"].encoding["dtype"] = "i4"
        ds_avg.to_netcdf(nc_filename, unlimited_dims=["time"])
        utils.check_compliance(nc_filename, conventions=ds_avg.attrs["Conventions"])

        print(f"Done writing netCDF file {nc_filename}")


def fill_time_gaps(ds):

    # Get time difference between each row
    time_diffs = np.diff(ds.time.values)

    # Check for more than one unique value
    all_same = np.unique(time_diffs).size == 1

    # If more than one unique value, fill gaps in time
    if not all_same:

        delta_t = int(time_diffs[0])

        # Make new timestamp
        timestamp = np.array(
            pd.date_range(
                start=ds.time[0].values, end=ds.time[-1].values, freq=f"{delta_t}ns"
            )
        )

        # Reindex onto continuous time
        ds = ds.reindex(time=timestamp, method="nearest")

    return ds


def reshape_var(ds, new_shape, new_dims):
    # Reshape values and put back into xarray
    reshaped_values = ds.values.reshape(new_shape)
    ds_reshape = xr.DataArray(reshaped_values, dims=new_dims)

    return ds_reshape


def reshape_burst(ds):
    # Reshape burst data because it is in one column

    # Calculate how many columns (samples) for each burst
    cols = int(ds.attrs["burst_duration"] / ds.attrs["sample_interval"])

    # If incomplete burst (time doesn't divide equally by samples_per_burst), remove remainder (trim incomplete burst)
    remainder = len(ds["time"]) % cols
    if remainder != 0:
        ds = ds.isel(time=slice(0, -remainder))

    # Calculate how many rows needed
    rows = int(len(ds["time"]) / cols)

    # Make new timestamp
    date = np.array(
        pd.date_range(
            start=ds.time[0].values, periods=rows, freq=f"{ds.attrs['burst_interval']}s"
        )
    )

    # Reshape
    new_shape = (rows, cols)
    new_dims = ("time", "sample")
    new_coords = {"time": date, "sample": np.arange(1, cols + 1)}

    ds_burst = xr.Dataset(
        {
            name: reshape_var(var, new_shape, new_dims)
            for name, var in ds.data_vars.items()
        },
        coords=new_coords,
    )

    ds_burst.attrs = ds.attrs

    return ds_burst


def avg_cont(ds):
    # Average continuous data

    # Calculate how many columns (samples) for each burst
    cols = int(ds.attrs["average_interval"] / ds.attrs["sample_interval"])

    # Calculate how many rows to subdivide continuous data into bursts
    rows = math.ceil(ds["time"].size / cols)

    # Make new timestamp
    date = np.array(
        pd.date_range(
            start=ds.time[0].values,
            periods=rows,
            freq=f"{ds.attrs['average_interval']}s",
        )
    )

    # Add padding to end of single column of continuous data if final burst will be incomplete and reshape into bursts
    no_pads = rows * cols - ds["time"].size

    temp_burst = np.pad(
        ds.T_28.values, (0, no_pads), mode="constant", constant_values=np.nan
    ).reshape(rows, cols)
    u_burst = np.pad(
        ds.u_1205.values, (0, no_pads), mode="constant", constant_values=np.nan
    ).reshape(rows, cols)
    v_burst = np.pad(
        ds.v_1206.values, (0, no_pads), mode="constant", constant_values=np.nan
    ).reshape(rows, cols)

    # Take mean of whole interval if no duration specified
    if "average_duration" not in ds.attrs:
        temp_avg = np.mean(temp_burst, axis=1)
        u_avg = np.mean(u_burst, axis=1)
        v_avg = np.mean(v_burst, axis=1)

        ds.attrs["samples_per_burst"] = cols

    # Take mean of duration
    elif "average_duration" in ds.attrs:

        # Calculate number of samples to average based on duration
        samp_per_burst = int(
            float(ds.attrs["average_duration"]) * float(ds.attrs["sample_rate"])
        )

        ds.attrs["samples_per_burst"] = samp_per_burst

        # Average burst for duration
        temp_avg = []
        u_avg = []
        v_avg = []
        for j in range(0, temp_burst.shape[0]):
            avg_temp_burst = np.mean(temp_burst[j, slice(0, samp_per_burst)])
            avg_u_burst = np.mean(u_burst[j, slice(0, samp_per_burst)])
            avg_v_burst = np.mean(v_burst[j, slice(0, samp_per_burst)])

            temp_avg.append(avg_temp_burst)
            u_avg.append(avg_u_burst)
            v_avg.append(avg_v_burst)

    # Make new xarray with averaged data
    ds_avg = xr.Dataset(
        data_vars=dict(
            T_28=(["time"], temp_avg),
            u_1205=(["time"], u_avg),
            v_1206=(["time"], v_avg),
        ),
        coords=dict(time=date),
    )

    return ds_avg


def ds_rename_vars(ds):
    """
    Rename variables to be CF compliant
    """
    varnames = {
        "temp": "T_28",
        "speed": "CS_300",
        "heading": "CD_310",
        "u": "u_1205",
        "v": "v_1206",
    }

    # Check to make sure they exist before trying to rename
    newvars = {}
    for k in varnames:
        if k in ds:
            newvars[k] = varnames[k]
    return ds.rename(newvars)


def ds_add_attrs(ds):
    """
    Add attributes: units, standard name from CF website, long names
    """

    ds["time"].attrs.update(
        {"standard_name": "time", "axis": "T", "long_name": "time (UTC)"},
    )

    if "T_28" in ds:
        ds["T_28"].attrs.update(
            {
                "units": "degree_C",
                "standard_name": "sea_water_temperature",
                "long_name": "Temperature",
            }
        )

    if "C_51" in ds:
        ds["C_51"].attrs.update(
            {
                "units": "S/m",
                "long_name": "Conductivity",
                "standard_name": "sea_water_electrical_conductivity",
            }
        )

    if "S_41" in ds:
        ds["S_41"].attrs.update(
            {
                "units": "1",
                "long_name": "Salinity, PSU",
                "comments": "Practical salinity units (PSU)",
                "standard_name": "sea_water_practical_salinity",
            }
        )

    if "u_1205" in ds:
        ds["u_1205"].attrs.update(
            {
                "units": "m s^-1",
                "long_name": "Eastward Velocity",
            }
        )

    if "v_1206" in ds:
        ds["v_1206"].attrs.update(
            {
                "units": "m s^-1",
                "long_name": "Northward Velocity",
            }
        )

    if "CS_300" in ds:
        ds["CS_300"].attrs.update(
            {
                "units": "m s^-1",
                "long_name": "Current Speed",
                "standard_name": "sea_water_speed",
            }
        )

    if "CD_310" in ds:
        ds["CD_310"].attrs.update(
            {
                "units": "degree",
                "long_name": "Current Direction (True)",
                "standard_name": "sea_water_velocity_to_direction",
            }
        )

    if "sample" in ds:
        ds["sample"].attrs.update(
            {
                "units": "1",
                "long_name": "sample number",
            }
        )
    return ds
