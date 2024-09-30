import pandas as pd
import xarray as xr

from .core import filter, qaqc, utils


def read_tid(filnam, encoding="utf-8"):
    """Read data from an SBE 26plus Seagauge .tid file into an xarray dataset."""
    df = pd.read_csv(
        filnam,
        header=None,
        sep=r"\s+",
        names=["Sample", "Date", "Time", "P_1", "Temp"],
        parse_dates={"time": ["Date", "Time"]},
        encoding=encoding,
        index_col=False,
    )
    df.set_index("time", inplace=True)
    sg = df.to_xarray()
    return sg


def tid_to_cdf(metadata):
    """
    Load a raw .tid and .hex file and generate a .cdf file
    """
    basefile = metadata["basefile"]

    # Get metadata from .hex file
    hexmeta = read_hex(basefile + ".hex")

    # Append to metadata variable
    metadata.update(hexmeta)

    # Read in data
    ds = read_tid(basefile + ".tid")

    # Convert pressure from psia to dbar
    ds["P_1"] = ds.P_1 / 14.503773800722 * 10

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

    # Drop sample variable
    ds = ds.drop_vars("Sample")

    # Rename variables to CF compliant names
    ds = ds_rename_vars(ds)

    # Atmospheric pressure correction
    if atmpres is not None:
        ds = atmos_correct(ds, atmpres)

    # Add attributes
    ds = ds_add_attrs(ds)

    # Call QAQC
    ds = sg_qaqc(ds)

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
    nc_filename = ds.attrs["filename"] + "-a.nc"

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
    return ds


def sg_qaqc(ds):
    """
    QA/QC
    Trim Seagauge data based on metadata
    """

    varlist = ["T_28", "P_1", "P_1ac"]

    [varlist.append(k) for k in ds.data_vars if k not in varlist]

    for var in varlist:
        # check if any filtering before other qaqc
        ds = filter.apply_butter_filt(ds, var)
        ds = filter.apply_med_filt(ds, var)

        ds = qaqc.trim_min(ds, var)

        ds = qaqc.trim_max(ds, var)

        ds = qaqc.trim_min_diff(ds, var)

        ds = qaqc.trim_min_diff_pct(ds, var)

        ds = qaqc.trim_max_diff(ds, var)

        ds = qaqc.trim_max_diff_pct(ds, var)

        ds = qaqc.trim_med_diff(ds, var)

        ds = qaqc.trim_med_diff_pct(ds, var)

        ds = qaqc.trim_bad_ens(ds, var)

    for var in varlist:
        ds = qaqc.trim_by_any(
            ds, var
        )  # re-run and trim by other variables as necessary

    return ds


def atmos_correct(ds, atmpres):
    met = xr.load_dataset(atmpres)
    # need to save attrs before the subtraction, otherwise they are lost
    attrs = ds["P_1"].attrs
    # need to set a tolerance since we can be off by a couple seconds somewhere
    # TODO is this still needed?
    ds["P_1ac"] = xr.DataArray(
        ds["P_1"]
        - met["atmpres"].reindex_like(ds["P_1"], method="nearest", tolerance="5s")
        - met["atmpres"].attrs["offset"]
    )
    ds["P_1ac"].attrs = attrs

    ds.attrs["atmospheric_pressure_correction_file"] = atmpres
    ds.attrs["atmospheric_pressure_correction_offset_applied"] = met["atmpres"].attrs[
        "offset"
    ]

    histtext = f"Atmospherically corrected using time-series from {atmpres} and offset of {met['atmpres'].offset}"

    ds = utils.insert_history(ds, histtext)

    # Also add it as a note
    ds["P_1ac"].attrs["note"] = histtext

    return ds


def read_hex(filnam):
    """Read metadata from an SBE 26plus Seagauge .hex file"""
    hexmeta = {}

    f = open(filnam)
    row = ""

    while "S>DD" not in row:
        row = f.readline().rstrip()
        if "Software Version" in row:
            col = row.split()
            hexmeta["SoftwareVersion"] = col[2]
        elif "SBE 26plus-quartz V" in row:
            col = row.split()
            hexmeta["InstrumentType"] = col[0][1:] + " " + col[1]
            hexmeta["FirmwareVersion"] = col[3]
            hexmeta["InstrumentSerialNumber"] = col[5]
        elif "quartz pressure sensor" in row:
            col = row.split()
            hexmeta["PressureSensorSerial"] = col[6][:-1]
        elif "tide measurement: interval" in row:
            col = row.split()
            hexmeta["TideInterval"] = col[4]
            hexmeta["TideIntervalUnits"] = col[5][:-1]
            hexmeta["TideDuration"] = col[8]
            hexmeta["TideDurationUnits"] = col[9]
        elif "measure waves every" in row:
            col = row.split()
            hexmeta["WaveInterval"] = col[3]
            hexmeta["WaveIntervalUnits"] = col[4] + " " + col[5]
        elif "wave samples/burst" in row:
            col = row.split()
            hexmeta["WaveSamples"] = col[0][1:]
            hexmeta["WaveSampleRate"] = col[4]
            hexmeta["WaveSampleRateUnits"] = col[5][:-1]
            hexmeta["BurstDuration"] = col[8]
            hexmeta["BurstDurationUnits"] = col[9]
        elif "tide samples/day" in row:
            col = row.split()
            hexmeta["TideSamplesPerDay"] = col[3]
        elif "wave bursts/day" in row:
            col = row.split()
            hexmeta["WaveBurstsPerDay"] = col[3]
        elif "total recorded tide measurements" in row:
            col = row.split()
            hexmeta["NumberOfTideMeasurements"] = col[5]
        elif "total recorded wave bursts" in row:
            col = row.split()
            hexmeta["NumberOfWaveBursts"] = col[5]
        elif "Pressure coefficients" in row:
            col = row.split()
            hexmeta["PressureCalibrationDate"] = col[2]
        elif "U0 =" in row:
            col = row.split()
            hexmeta["PressureCalibrationU0"] = float(col[3])
        elif "Y1 =" in row:
            col = row.split()
            hexmeta["PressureCalibrationY1"] = float(col[3])
        elif "Y2 =" in row:
            col = row.split()
            hexmeta["PressureCalibrationY2"] = float(col[3])
        elif "Y3 =" in row:
            col = row.split()
            hexmeta["PressureCalibrationY3"] = float(col[3])
        elif "C1 =" in row:
            col = row.split()
            hexmeta["PressureCalibrationC1"] = float(col[3])
        elif "C2 =" in row:
            col = row.split()
            hexmeta["PressureCalibrationC2"] = float(col[3])
        elif "C3 =" in row:
            col = row.split()
            hexmeta["PressureCalibrationC3"] = float(col[3])
        elif "D1 =" in row:
            col = row.split()
            hexmeta["PressureCalibrationD1"] = float(col[3])
        elif "D2 =" in row:
            col = row.split()
            hexmeta["PressureCalibrationD2"] = float(col[3])
        elif "T1 =" in row:
            col = row.split()
            hexmeta["PressureCalibrationT1"] = float(col[3])
        elif "T2 =" in row:
            col = row.split()
            hexmeta["PressureCalibrationT2"] = float(col[3])
        elif "T3 =" in row:
            col = row.split()
            hexmeta["PressureCalibrationT3"] = float(col[3])
        elif "T4 =" in row:
            col = row.split()
            hexmeta["PressureCalibrationT4"] = float(col[3])
        elif "M =" in row:
            col = row.split()
            hexmeta["PressureCalibrationM"] = float(col[3])
        elif "B =" in row:
            col = row.split()
            hexmeta["PressureCalibrationB"] = float(col[3])
        elif "OFFSET =" in row:
            col = row.split()
            hexmeta["PressureCalibrationOFFSET"] = float(col[3])
        elif "Temperature coefficients" in row:
            col = row.split()
            hexmeta["TemperatureCalibrationDate"] = col[2]
        elif "TA0 =" in row:
            col = row.split()
            hexmeta["TemperatureCalibrationTA0"] = float(col[3])
        elif "TA1 =" in row:
            col = row.split()
            hexmeta["TemperatureCalibrationTA1"] = float(col[3])
        elif "TA2 =" in row:
            col = row.split()
            hexmeta["TemperatureCalibrationTA2"] = float(col[3])
        elif "TA3 =" in row:
            col = row.split()
            hexmeta["TemperatureCalibrationTA3"] = float(col[3])
    f.close()
    return hexmeta
