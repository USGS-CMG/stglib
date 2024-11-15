from datetime import timedelta

import pandas as pd
import xarray as xr

from .core import filter, qaqc, utils


def read_header(filnam):
    """Read header info from a TruBlue 255 .csv file to add to metadata."""

    header = {}

    f = open(filnam)
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

    f.close()
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
    df["time"] = df["time"].astype("datetime64[ns]")
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


def cdf_to_nc(cdf_filename, atmpres=None):
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
    ds = ds.drop("ID")
    ds = ds.drop("Name")
    ds = ds.drop("Address")
    ds = ds.drop("Elapsed")

    if atmpres:
        ds = utils.atmos_correct(ds, atmpres)

    # Add attributes
    ds = ds_add_attrs(ds)

    # Call QAQC
    ds = tb_qaqc(ds)

    # Run utilities
    ds = utils.clip_ds(ds)
    ds = utils.create_nominal_instrument_depth(ds)
    ds = utils.create_z(ds)
    ds = utils.create_water_level_var(ds)
    ds = utils.create_water_depth_var(ds)
    ds = utils.ds_add_lat_lon(ds)
    ds = utils.add_start_stop_time(ds)
    ds = utils.add_delta_t(ds)
    ds = utils.add_min_max(ds)

    # Write to .nc file
    print("Writing cleaned/trimmed data to .nc file")
    nc_filename = ds.attrs["filename"] + "-cont-cal.nc"

    ds.to_netcdf(
        nc_filename, unlimited_dims=["time"], encoding={"time": {"dtype": "i4"}}
    )
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])

    print(f"Done writing netCDF file {nc_filename}")


def ds_add_attrs(ds):
    """
    Add attributes: units, standard name from CF website, long names
    """
    ds = utils.ds_coord_no_fillvalue(ds)

    ds["time"].attrs.update(
        {"standard_name": "time", "axis": "T", "long_name": "time (UTC)"}
    )

    if "T_28" in ds:
        ds["T_28"].attrs.update(
            {
                "units": "degree_C",
                "standard_name": "sea_water_temperature",
                "long_name": "Temperature",
            }
        )

    if "P_1" in ds:
        ds["P_1"].attrs.update(
            {
                "units": "dbar",
                "long_name": "Uncorrected pressure",
                "standard_name": "sea_water_pressure",
            }
        )

    if "P_1ac" in ds:
        ds["P_1ac"].attrs.update(
            {
                "units": "dbar",
                "long_name": "Corrected pressure",
                "standard_name": "sea_water_pressure_due_to_sea_water",
            }
        )
        if "P_1ac_note" in ds.attrs:
            ds["P_1ac"].attrs.update({"note": ds.attrs["P_1ac"]})

    return ds


def tb_qaqc(ds):
    """
    QA/QC
    Trim data based on metadata
    """

    varlist = ["T_28", "P_1", "P_1ac"]

    for var in varlist:
        if var in ds:
            ds = filter.apply_butter_filt(ds, var)
            ds = filter.apply_med_filt(ds, var)

            ds = qaqc.trim_bad_ens(ds, var)
            ds = qaqc.trim_bad_ens_indiv(ds, var)

            ds = qaqc.trim_min(ds, var)
            ds = qaqc.trim_max(ds, var)

            ds = qaqc.trim_min_diff(ds, var)
            ds = qaqc.trim_min_diff_pct(ds, var)

            ds = qaqc.trim_max_diff(ds, var)
            ds = qaqc.trim_max_diff_pct(ds, var)

            ds = qaqc.trim_med_diff(ds, var)
            ds = qaqc.trim_med_diff_pct(ds, var)

            ds = qaqc.trim_max_blip(ds, var)
            ds = qaqc.trim_max_blip_pct(ds, var)

            ds = qaqc.trim_fliers(ds, var)
            ds = qaqc.trim_mask(ds, var)

            ds = qaqc.trim_maxabs_diff(ds, var)
            ds = qaqc.drop_vars(ds)

    for var in varlist:
        if var in ds:
            ds = qaqc.trim_by_any(
                ds, var
            )  # re-run and trim by other variables as necessary

    return ds
