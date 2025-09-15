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

    # Log file with metadata. Only one log file even if multiple exported ASCII and .mat files
    ds = read_log_file(ds, basefile.with_suffix(".log"))

    # Matlab file with exported velocity data (but no pressure data, which is why we need the above).
    matfiles = [p for p in basefile.parent.iterdir() if p.suffix == ".mat"]
    # mat = read_mat_file(basefile.with_suffix(".000.mat"))
    mats = read_mat_file(matfiles)

    # File with exported pressure data
    # sens = read_sens_file(basefile.with_suffix(".000.txt"))
    sens = read_sens_file([Path(str(p).replace(".mat", ".txt")) for p in matfiles])

    dss = []
    for mat in mats:
        dsm = xr.Dataset()
        dsm["time"] = pd.to_datetime(mat["sens"]["time"], unit="s")
        dsm["time"].attrs["standard_name"] = "time"
        dsm["time"].encoding["dtype"] = "int32"

        dsm["bindist"] = mat["info"]["cell1"] + mat["info"]["cell"] * np.arange(
            int(mat["info"]["ncells"])
        )
        dsm["velbeam"] = ["E", "N", "U1", "U2"]
        dsm["velbeam"].encoding["dtype"] = "str"
        dsm["beam"] = [1, 2, 3, 4]
        dsm["beam"].encoding["dtype"] = "i4"

        dsm["bindist"].attrs["units"] = "m"
        dsm["bindist"].attrs[
            "long_name"
        ] = "bin distance from instrument for slant beams"
        dsm["bindist"].attrs["epic_code"] = 0
        dsm["bindist"].attrs[
            "NOTE"
        ] = "distance is calculated from center of bin 1 and bin size"
        dsm["vel"] = (["time", "bindist", "velbeam"], mat["wt"]["vel"])

        for k in ["int", "corr", "pg"]:
            dsm[k] = (["time", "bindist", "beam"], mat["wt"][k])

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
            dsm[sensnames[k]] = ("time", mat["sens"][k])

        for k in mat["info"]:
            dsm.attrs[k] = mat["info"][k]

        dsm["Heading"].attrs["units"] = "degrees"
        dsm["Heading"].attrs["standard_name"] = "platform_orientation"
        dsm["Heading"].attrs["epic_code"] = 1215

        dsm["Pitch"].attrs["units"] = "degrees"
        dsm["Pitch"].attrs["standard_name"] = "platform_pitch"
        dsm["Pitch"].attrs["epic_code"] = 1216

        dsm["Roll"].attrs["units"] = "degrees"
        dsm["Roll"].attrs["standard_name"] = "platform_roll"
        dsm["Roll"].attrs["epic_code"] = 1217

        dsm["sv"].attrs["units"] = "m s-1"
        dsm["sv"].attrs["standard_name"] = "speed_of_sound_in_sea_water"

        dsm["S"].attrs["units"] = "1"
        dsm["S"].attrs["standard_name"] = "sea_water_salinity"
        dsm["S"].attrs["epic_code"] = 40

        dsm["Temperature"].attrs["units"] = "degree_C"
        dsm["Temperature"].attrs["long_name"] = "ADCP Transducer Temperature"
        dsm["Temperature"].attrs["epic_code"] = 1211

        if "vel5" in ds:  # 5-beam instrument
            dsm["vel5"].attrs["units"] = "mm s-1"
            dsm["vel5"].attrs["long_name"] = "Beam 5 velocity (mm s-1)"
            dsm["vel5"].encoding["dtype"] = "i2"

            dsm["cor5"].attrs["units"] = "counts"
            dsm["cor5"].attrs["long_name"] = "Beam 5 correlation"
            dsm["cor5"].encoding["dtype"] = "u2"

            dsm["att5"].attrs["units"] = "counts"
            dsm["att5"].attrs["long_name"] = "ADCP attenuation of beam 5"
            dsm["att5"].encoding["dtype"] = "u2"
        dss.append(dsm)

    dsm = xr.concat(dss, dim="time", combine_attrs="no_conflicts")

    ds = xr.merge([ds, dsm], combine_attrs="no_conflicts")

    if "Pressure" in sens:
        ds["Pressure"] = sens["Pressure"] / 10  # convert to dbar
        ds["Pressure"].attrs["units"] = "dbar"
        ds["Pressure"].attrs["standard_name"] = "sea_water_pressure"

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
    if not isinstance(filnam, list):
        filnam = [filnam]

    mats = []
    for f in filnam:
        mats.append(utils.loadmat(f))
    return mats


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
    if not isinstance(filnam, list):
        filnam = [filnam]
    senss = []
    for f in filnam:
        sens = pd.read_csv(f, header=0)

        # Need to be named "minute", "second" for to_datetime() to work
        sens = sens.rename(columns={"Min": "Minute", "Sec": "Second"})
        sens["time"] = pd.to_datetime(
            sens[["Year", "Month", "Day", "Hour", "Minute", "Second"]]
        )
        senss.append(sens)

    return pd.concat(senss).set_index("time")
