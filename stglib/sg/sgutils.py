import math

import numpy as np
import pandas as pd
import xarray as xr

from ..core import qaqc, utils


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
            hexmeta["serial_number"] = col[5]
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
    start_time = int_to_date(np.array(start_time))

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


def int_to_date(int_time):
    """Read integer time in seconds since January 1, 2000 and convert to datetime"""

    dt = pd.Timestamp("2000-01-01T00:00:00") + pd.to_timedelta(int_time, unit="s")

    # This gave me a nanosecond precision error
    # t0 = np.datetime64("2000-01-01T00:00:00")
    # dt = t0 + int_time.astype("timedelta64[s]")

    return dt


def atmos_correct_burst(ds, atmpres):
    met = xr.load_dataset(atmpres)
    pressure = []

    # Apply the correction for each burst in turn
    ds["P_1ac"] = xr.full_like(ds["P_1"], np.nan)
    for burst in ds.burst_number:
        burst_pres = (
            ds["P_1"][burst].values
            - met["atmpres"][burst].values
            - met["atmpres"].offset
        )
        burst_pres = np.reshape(burst_pres, (1, -1))

        # Convert back to list to append bursts
        burst_pres = burst_pres.tolist()
        pressure.append(burst_pres)

    # Convert to xarray
    pressure = xr.DataArray.squeeze(
        xr.DataArray(pressure, dims=["time", "one", "sample"], name="P_1ac")
    )
    ds = xr.merge([ds, pressure])

    ds = utils.insert_history(
        ds,
        f"Atmospherically correcting using time-series from {atmpres} and offset of {met['atmpres'].offset}",
    )
    ds.attrs["atmospheric_pressure_correction_file"] = atmpres
    ds.attrs["atmospheric_pressure_correction_offset_applied"] = met["atmpres"].attrs[
        "offset"
    ]
    if "comment" in met["atmpres"].attrs:
        ds.attrs["atmospheric_pressure_correction_comment"] = met["atmpres"].attrs[
            "comment"
        ]

    return ds


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
    if "sample" in ds:
        ds["sample"].attrs.update(
            {
                "units": "1",
                "long_name": "sample number",
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
