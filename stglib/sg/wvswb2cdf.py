import math

import numpy as np
import xarray as xr

from ..core import utils
from . import sgutils


def read_wb(filnam, encoding="utf-8"):
    """Read data from an SBE 26plus Seagauge .wb file into an xarray dataset."""

    burst_no = []
    start_time = []
    values = []
    pressure = []

    with open(filnam) as file:
        for line in file:
            if "SBE" in line:
                continue
            elif "*" in line:
                col = line.split()
                burst_no.append(int(col[1]))
                start_time.append(int(col[2]))
                sample_no = int(col[4])
                sample_rows = math.floor(
                    sample_no / 4
                )  # Take total number of samples in each burst and divide by 4 columns to get number of rows to iterate

                # Loop through rows for each burst number
                for j in range(sample_rows):

                    row = file.readline().rstrip()
                    col = row.split()
                    values.append(col)

                # Convert to numpy float and reshape to 1 row
                burst = np.float64(values)
                burst = np.reshape(burst, (1, -1))

                # Convert back to list to append bursts
                burst = burst.tolist()
                pressure.append(burst)

                # Clear values for next burst
                values = []
    file.close()

    # Convert pressure to numpy array
    pressure = np.float64(pressure)
    pressure = pressure.squeeze()

    # Convert time
    start_time = int_to_date(start_time)

    # Create sample variable
    sample = list(range(1, sample_no + 1))

    # Make xarray
    sg = xr.Dataset(
        data_vars=dict(
            burst_number=(["time"], burst_no),
            sample=(["sample"], sample),
            P_1=(["time", "sample"], pressure),
        ),
        coords=dict(time=start_time),
    )

    return sg

    # sg=xr.DataArray(pressure, coords=[start_time, sample], dims=['time','sample'])

    # sg=xr.DataArray({"Pressure":pressure, "time":start_time, "sample": sample}, coords=[start_time, sample], dims=['time','sample'])

    # pr["row_start"] = xr.DataArray(row_start, dims="profile")
    # xr.merge


def int_to_date(int_time):
    """Read integer time in seconds since January 1, 2000 and convert to datetime"""

    t0 = np.datetime64("2000-01-01T00:00:00")

    dt = t0 + int_time.astype("timedelta64[s]")

    return dt


def wb_to_cdf(metadata):
    """
    Load a raw .wb and .hex file and generate a .cdf file
    """
    basefile = metadata["basefile"]

    # Get metadata from .hex file
    hexmeta = sgutils.read_hex(basefile + ".hex")

    # Append to metadata variable
    metadata.update(hexmeta)

    # Read in data
    ds = read_wb(basefile + ".wb")

    # Convert pressure from psia to dbar
    ds["P_1"] = ds.P_1 / 14.503773800722 * 10

    ds = utils.write_metadata(ds, metadata)

    ds = utils.ensure_cf(ds)

    cdf_filename = ds.attrs["filename"] + "-waves-raw.cdf"

    ds.to_netcdf(cdf_filename, unlimited_dims=["time"])

    print(f"Finished writing data to {cdf_filename}")

    return ds
