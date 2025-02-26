from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from ..aqd import aqdutils
from ..core import utils


def mat_to_cdf(metadata):
    """Load RDI ADCP text files and output to netCDF format"""

    utils.check_valid_globalatts_metadata(metadata)
    aqdutils.check_valid_config_metadata(metadata)

    ds = xr.Dataset()

    # write out metadata first, then deal exclusively with xarray attrs
    ds = utils.write_metadata(ds, metadata)

    basefile = Path(ds.attrs["basefile"])

    if "prefix" in ds.attrs:
        basefile = ds.attrs["prefix"] / basefile

    # Log file with metadata
    ds = read_log_file(ds, basefile.with_suffix(".log"))

    # File with exported pressure data
    sens = read_sens_file(basefile.with_suffix(".000.txt"))

    # Matlab file with exported velocity data (but no pressure data, which is why we need the above).
    mat = read_mat_file(basefile.with_suffix(".000.mat"))

    ds["time"] = pd.to_datetime(mat["sens"]["time"], unit="s")
    ds["time"].attrs["standard_name"] = "time"
    ds["time"].encoding["dtype"] = "int32"

    ds["bindist"] = mat["info"]["cell1"] + mat["info"]["cell"] * np.arange(
        int(mat["info"]["ncells"])
    )
    ds["velbeam"] = ["E", "N", "U1", "U2"]
    ds["beam"] = [1, 2, 3, 4]
    ds["beam"].encoding["dtype"] = "i4"

    ds["bindist"].attrs["units"] = "m"
    ds["bindist"].attrs["long_name"] = "bin distance from instrument for slant beams"
    ds["bindist"].attrs["epic_code"] = 0
    ds["bindist"].attrs[
        "NOTE"
    ] = "distance is calculated from center of bin 1 and bin size"

    ds["vel"] = (["time", "bindist", "velbeam"], mat["wt"]["vel"])

    for k in ["int", "corr", "pg"]:
        ds[k] = (["time", "bindist", "beam"], mat["wt"][k])

    sensnames = {
        "h": "Heading",
        "p": "Pitch",
        "r": "Roll",
        "t": "Temperature",
        "sos": "sv",
        "s": "S",
        "o": "Orient",
        "v": "Battery",
    }
    for k in sensnames:
        ds[sensnames[k]] = ("time", mat["sens"][k])

    for k in mat["info"]:
        ds.attrs[k] = mat["info"][k]

    ds["Heading"].attrs["units"] = "degrees"
    ds["Heading"].attrs["standard_name"] = "platform_orientation"
    ds["Heading"].attrs["epic_code"] = 1215

    ds["Pitch"].attrs["units"] = "degrees"
    ds["Pitch"].attrs["standard_name"] = "platform_pitch"
    ds["Pitch"].attrs["epic_code"] = 1216

    ds["Roll"].attrs["units"] = "degrees"
    ds["Roll"].attrs["standard_name"] = "platform_roll"
    ds["Roll"].attrs["epic_code"] = 1217

    ds["sv"].attrs["units"] = "m s-1"
    ds["sv"].attrs["standard_name"] = "speed_of_sound_in_sea_water"

    ds["S"].attrs["units"] = "1"
    ds["S"].attrs["standard_name"] = "sea_water_salinity"
    ds["S"].attrs["epic_code"] = 40

    ds["Temperature"].attrs["units"] = "degree_C"
    ds["Temperature"].attrs["long_name"] = "ADCP Transducer Temperature"
    ds["Temperature"].attrs["epic_code"] = 1211

    if "Pressure" in sens:
        ds["Pressure"] = sens["Pressure"] / 10  # convert to dbar
        ds["Pressure"].attrs["units"] = "dbar"
        ds["Pressure"].attrs["standard_name"] = "sea_water_pressure"

    if "vel5" in ds:  # 5-beam instrument
        ds["vel5"].attrs["units"] = "mm s-1"
        ds["vel5"].attrs["long_name"] = "Beam 5 velocity (mm s-1)"
        ds["vel5"].encoding["dtype"] = "i2"

        ds["cor5"].attrs["units"] = "counts"
        ds["cor5"].attrs["long_name"] = "Beam 5 correlation"
        ds["cor5"].encoding["dtype"] = "u2"

        ds["att5"].attrs["units"] = "counts"
        ds["att5"].attrs["long_name"] = "ADCP attenuation of beam 5"
        ds["att5"].encoding["dtype"] = "u2"

    cdf_filename = Path(ds.attrs["filename"] + "-raw.cdf")

    if "prefix" in ds.attrs:
        cdf_filename = ds.attrs["prefix"] / cdf_filename

    ds.to_netcdf(cdf_filename, unlimited_dims=["time"])

    print("Finished writing data to %s" % cdf_filename)

    return ds


def read_mat_file(filnam):
    # Format of .mat file
    # array(['info (System/Setup Info):          ',
    #        '   cell - depth cell size [m]      ',
    #        '   blank - blank [m]               ',
    #        '   cell1 - first cell range [m]    ',
    #        '   ncells - number of cells        ',
    #        '   angle - beam angle [deg]        ',
    #        '                                   ',
    #        'sens (Sensors):                    ',
    #        '   dnum - time int Datenum format  ',
    #        '   time - seconds from 1970/01/01  ',
    #        '   h - heading [deg]               ',
    #        '   p - pitch [deg]                 ',
    #        '   r - roll [deg]                  ',
    #        '   t - temperature [C]             ',
    #        '   pd - pressure sensor depth [m]  ',
    #        '   sos - speed of sound [m/s]      ',
    #        '   s - salinity                    ',
    #        '   o - orientation                 ',
    #        '   v - voltage [V]                 ',
    #        '                                   ',
    #        'wt (Water profile):                ',
    #        '   vel - velocity [m/s]            ',
    #        '   int - intensity [counts]        ',
    #        '   corr - correlation [counts]     ',
    #        '   pg - percent good [%]           ',
    #        '   d - cells depths [m]            ',
    #        '   r - cells ranges [m]            '], dtype=object)

    return utils.loadmat(filnam)


def read_log_file(ds, filnam):
    with open(filnam) as f:
        for line in f:
            for m in [
                "File size",
                "Valid data",
                "Invalid data",
                "Record size",
                "First record number",
                "First record time",
                "Last record number",
                "Last record time",
                "Total records",
                "Missing records",
                "Bad BIT records",
                "Software version",
                "Firmware version",
                "System type",
                "Serial number",
                "Frequency",
                "Number of cells",
                "Cell size",
                "Blank",
                "Water mode",
                "Water pings",
            ]:
                if m in line[:20]:
                    attrname = "".join(m.title().split())
                    ds.attrs[f"RDI{attrname}"] = line.split("\t")[1].strip()

    return ds


def read_sens_file(filnam):
    sens = pd.read_csv(filnam, header=0)

    # Need to be named "minute", "second" for to_datetime() to work
    sens = sens.rename(columns={"Min": "Minute", "Sec": "Second"})
    sens["time"] = pd.to_datetime(
        sens[["Year", "Month", "Day", "Hour", "Minute", "Second"]]
    )

    return sens.set_index("time")
