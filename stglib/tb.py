from datetime import timedelta

import pandas as pd
import xarray as xr

from .core import attrs, qaqc, utils


def read_header(filnam):
    """Read header info from a TruBlue 255 .csv file to add to metadata."""

    header = {}

    with open(filnam) as f:
        row = ""
        while "ID,Name,Address" not in row:
            row = f.readline().rstrip()
            col = row.split()
            if "Transducer Model" in row:
                header["TransducerModel"] = col[3]
            elif "Transducer Serial" in row:
                header["serial_number"] = col[3]
            elif "Measure Parameters" in row:
                header["MeasureParameters"] = col[2]
            elif "Test Code" in row:
                header["TestCode"] = col[2]
            elif "Scan Type" in row:
                header["ScanType"] = col[2]
            elif "Periods" in row:
                header["Periods"] = col[1]
            elif "Memory Wrap" in row:
                header["MemoryWrap"] = col[2]
            elif "Broadcast Normal" in row:
                header["BroadcastNormalLog"] = col[3]
            elif "Broadcast Alarm" in row:
                header["BroadcastAlarmLog"] = col[3]
            elif "Alarm Scan" in row:
                header["AlarmScan"] = col[2]
            elif "Scanned Upto" in row:
                header["ScannedUpto"] = col[3]
            elif "Firmware Version" in row:
                header["FirmwareVersion"] = col[2]
            elif "TruWare" in row:
                header["TruWareVersion"] = col[1]

    return header


def read_csv(filnam, skiprows=20, encoding="utf-8"):
    """Read data from a TruBlue 255 .csv file into an xarray
    Dataset.
    """

    df = pd.read_csv(
        filnam,
        skiprows=skiprows,
        header=None,
        names=["ID", "Name", "Address", "time", "Elapsed", "P_1", "T_28"],
        encoding=encoding,
        index_col=False,
    )
    df["time"] = df["time"].str.strip("'")
    df["time"] = pd.to_datetime(df["time"])
    df.set_index("time", inplace=True)
    tb = df.to_xarray()
    return tb


def txt_to_cdf(metadata):
    """
    Load a raw .txt file and generate a .cdf file
    """
    basefile = metadata["basefile"]

    # Read header info from file
    header = read_header(basefile + ".csv")

    # Add attribute sample_mode continuous to use RBR nc2waves.py
    header["sample_mode"] = "CONTINUOUS"

    # Append to metadata variable
    metadata.update(header)

    # Read raw data
    ds = read_csv(basefile + ".csv", skiprows=metadata["skiprows"])

    # Convert pressure from psi to dbar
    ds["P_1"] = ds.P_1 / 14.503773800722 * 10

    metadata.pop("skiprows")

    ds = utils.write_metadata(ds, metadata)

    ds = utils.ensure_cf(ds)

    cdf_filename = ds.attrs["filename"] + "-raw.cdf"

    ds.to_netcdf(cdf_filename, unlimited_dims=["time"])

    print(f"Finished writing data to {cdf_filename}")

    return ds


def cdf_to_nc(cdf_filename, atmpres=None, salwtemp=None):
    """
    Load a raw .cdf file and generate a processed .nc file
    """

    # Load raw .cdf data
    ds = xr.open_dataset(cdf_filename)

    # remove units in case we change and we can use larger time steps
    ds.time.encoding.pop("units")

    # Add sample_interval to metadata
    tdelta = timedelta(
        hours=int(ds.attrs["Periods"][3:5]),
        minutes=int(ds.attrs["Periods"][6:8]),
        seconds=float(ds.attrs["Periods"][9:]),
    )
    ds.attrs["sample_interval"] = tdelta.total_seconds()

    # Drop unneeded variables
    droplist = ["ID", "Name", "Address", "Elapsed"]
    ds = ds.drop_vars(droplist)

    if atmpres:
        ds = utils.atmos_correct(ds, atmpres)

    # Add attributes
    ds = attrs.ds_add_attrs(ds)

    # Run utilities
    ds = utils.clip_ds(ds)
    ds = utils.create_nominal_instrument_depth(ds)
    ds = utils.create_z(ds)
    ds = utils.create_water_level_var(ds, salwtemp=salwtemp)
    ds = utils.create_filtered_water_level_var(ds)
    ds = utils.create_water_depth_var(ds, salwtemp=salwtemp)
    ds = qaqc.call_qaqc(ds)  # Needs to happen after water level
    ds = utils.ds_add_lat_lon(ds)
    ds = utils.add_start_stop_time(ds)
    ds = utils.add_min_max(ds)

    if ds.attrs["sample_interval"] >= 1:
        ds = utils.add_delta_t(ds)

    # Write to .nc file
    print("Writing cleaned/trimmed data to .nc file")
    nc_filename = ds.attrs["filename"] + "-cont-cal.nc"

    ds.to_netcdf(
        nc_filename, unlimited_dims=["time"], encoding={"time": {"dtype": "i4"}}
    )
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])

    print(f"Done writing netCDF file {nc_filename}")
