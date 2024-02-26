import csv
import string
import warnings

import numpy as np
import pandas as pd
import xarray as xr

from .aqd import aqdutils
from .core import qaqc, utils


def read_tcm(
    filnam,
    skiprows=1,
    skipfooter=0,
    names=["DateTime", "Speed", "Bearing", "Velocity-N", "Velocity-E"],
):
    """Read data from an Lowell Instruments Tilt Current Meter (TCM) _CR.txt file into an xarray
    Dataset.

    Parameters
    ----------
    filnam : string
        The filename
    skiprows : int, optional
        How many header rows to skip. Default 1
    skipfooter : int, optional
        How many footer rows to skip. Default 0

    Returns
    -------
    xarray.Dataset
        An xarray Dataset of the TCM data
    """
    tcm = pd.read_csv(
        filnam,
        usecols=np.arange(len(names)),
        names=names,
        engine="python",
        skiprows=skiprows,
        skipfooter=skipfooter,
    )

    tcm["time"] = pd.to_datetime(tcm["DateTime"])

    tcm.set_index("time", inplace=True)

    return xr.Dataset(tcm)


def csv_to_cdf(metadata):
    """
    Process TCM .txt file to a raw .cdf file
    """

    basefile = metadata["basefile"]

    kwargs = {"skiprows": metadata["skiprows"], "skipfooter": metadata["skipfooter"]}
    if "names" in metadata:
        kwargs["names"] = metadata["names"]
    else:
        kwargs["names"] = get_tcm_col_names(basefile + "_CR.txt", metadata)

    try:
        ds = read_tcm(basefile + "_CR.txt", **kwargs)
    except IOError:
        print(f"Could not read file {basefile}_CR.txt, check file encoding, use utf-8")

    metadata.pop("skiprows")
    metadata.pop("skipfooter")
    if "ncols" in metadata:
        metadata.pop("ncols")

    # write out metadata first, then deal exclusively with xarray attrs
    ds = utils.write_metadata(ds, metadata)

    del metadata

    ds = utils.ensure_cf(ds)

    ds = utils.shift_time(ds, 0)

    ds = drop_vars(ds)

    ds.attrs["serial_number"] = get_serial_number(basefile + ".lid")

    # configure file
    cdf_filename = ds.attrs["filename"] + "-raw.cdf"

    ds.to_netcdf(cdf_filename, unlimited_dims=["time"])

    print("Finished writing data to %s" % cdf_filename)

    return ds


def drop_vars(ds):
    todrop = ["DateTime"]
    return ds.drop([x for x in todrop if x in ds])


def ds_rename_vars(ds):
    # convert some units

    if "Velocity-E_cms" in ds:
        ds["Velocity-E_cms"].values = (
            ds["Velocity-E_cms"].values / 100
        )  # convert from cm/s to m/s
        ds = ds.rename({"Velocity-E_cms": "u_1205"})

    if "Velocity-N_cms" in ds:
        ds["Velocity-N_cms"].values = (
            ds["Velocity-N_cms"].values / 100
        )  # convert from cm/s to m/s
        ds = ds.rename({"Velocity-N_cms": "v_1206"})

    if "Speed_cms" in ds:
        ds["Speed_cms"].values = (
            ds["Speed_cms"].values / 100
        )  # convert from cm/s to m/s
        ds = ds.rename({"Speed_cms": "CS_300"})

    if "Bearing_degrees" in ds:
        ds = ds.rename({"Bearing_degrees": "CD_310"})

    # drop unneeded vars
    todrop = [""]
    ds = ds.drop([x for x in todrop if x in ds])

    return ds


def ds_add_attrs(ds):
    # Update attributes for EPIC and STG compliance
    ds = utils.ds_coord_no_fillvalue(ds)

    ds["time"].attrs.update(
        {"standard_name": "time", "axis": "T", "long_name": "time (UTC)"}
    )

    if "u_1205" in ds:
        ds["u_1205"].attrs.update(
            {
                "units": "m s^-1",
                "long_name": "Eastward Velocity",
                "epic_code": 1205,
            }
        )

    if "v_1206" in ds:
        ds["v_1206"].attrs.update(
            {
                "units": "m s^-1",
                "long_name": "Northward Velocity",
                "epic_code": 1206,
            }
        )

    if "T_28" in ds:
        ds["T_28"].attrs.update(
            {
                "units": "degree_C",
                "long_name": "Temperature",
                "epic_code": 28,
                "standard_name": "sea_water_temperature",
            }
        )

    if "CS_300" in ds:
        ds["CS_300"].attrs.update(
            {
                "units": "m s^-1",
                "long_name": "Current Speed",
                "epic_code": 300,
                "standard_name": "sea_water_speed",
            }
        )

    if "CD_310" in ds:
        ds["CD_310"].attrs.update(
            {
                "units": "degree",
                "long_name": "Current Direction (True)",
                "epic_code": 310,
                "standard_name": "sea_water_velocity_to_direction",
            }
        )

    return ds


def get_serial_number(filnam):
    """get the serial number of the instrument from the lid file"""
    with open(filnam, "rb") as f:
        while True:
            row = f.readline()
            print(row)
            if b"SER " in row:
                sn = row.split()[1]
                return sn.decode("ascii")


def strip_non_printable(strin):
    """Returns the string without non printable characters"""
    printable = set(string.printable)
    return "".join(filter(lambda x: x in printable, strin))


def get_tcm_col_names(filnam, metadata):
    """get column names and column units from instrument input data file"""

    with open(filnam) as f:
        rdr = csv.reader(f)
        # check to see if "Time" in first
        hdrline = next(rdr)
        while "Time" not in hdrline[0]:
            hdrline = next(rdr)

    collist = [x.split(")")[0].replace(" ", "") for x in hdrline]
    # set first column to DateTime
    collist[0] = "DateTime"
    # print(collist)

    colnames = []
    colunits = []
    for x in collist:
        # spl = x.split(",")
        spl = x.replace(" ", "").split("(")
        colnames.append(spl[0].strip("."))
        if len(spl) > 1:
            colunits.append(spl[1].strip())
        else:
            colunits.append("")

    if "ncols" in metadata:
        colnames = colnames[: metadata["ncols"]]
        colunits = colunits[: metadata["ncols"]]

    # make dict of names and units
    dcols = {}
    for i in range(len(colnames)):
        d = {colnames[i]: colunits[i]}
        dcols.update(d)

    # try removing special characters and those not allowed in var or dim names from units
    for k in dcols:
        # first step try replacing values
        if "µ" in dcols[k]:
            dcols[k] = dcols[k].replace("µ", "u")
        if "°" in dcols[k]:
            dcols[k] = dcols[k].replace("°", "")
        if "%" in dcols[k]:
            dcols[k] = dcols[k].replace("%", "percent")
        if "Temp" in k:
            if "C" in dcols[k]:
                dcols[k] = "C"
            elif "F" in dcols[k]:
                dcols[k] = "F"
        if "/" in dcols[k]:
            dcols[k] = dcols[k].replace("/", "")

        # then strip non-ascii characters
        dcols[k] = strip_non_printable(dcols[k])

    names = []
    for k in dcols:
        if k == "#" or k == "DateTime":
            names.append(k)
        else:
            names.append(k + "_" + dcols[k])

    return names


def magvar_correct(ds):
    """Correct for magnetic declination at site"""

    if "magnetic_variation_at_site" in ds.attrs:
        magvardeg = ds.attrs["magnetic_variation_at_site"]
    elif "magnetic_variation" in ds.attrs:
        magvardeg = ds.attrs["magnetic_variation"]
    else:
        warnings.warn(
            "No magnetic variation information provided; not correcting compass orientation"
        )
        magvardeg = 0

    if magvardeg != 0:
        histtext = "Rotating heading and horizontal velocities by {} degrees.".format(
            magvardeg
        )

        ds = utils.insert_history(ds, histtext)

        if "U" in ds and "V" in ds:
            uvar = "U"
            vvar = "V"
        elif "u_1205" in ds and "v_1206" in ds:
            uvar = "u_1205"
            vvar = "v_1206"

        ds[uvar], ds[vvar] = aqdutils.rotate(ds[uvar], ds[vvar], magvardeg)

        if "CD_310" in ds:
            dvar = "CD_310"

        ds[dvar] = (ds[dvar] + magvardeg) % 360

    return ds


def cdf_to_nc(cdf_filename):
    """
    Load a raw .cdf file and generate a processed .nc file
    """

    # Load raw .cdf data
    ds = xr.open_dataset(cdf_filename)

    # remove units in case we change and we can use larger time steps
    ds.time.encoding.pop("units")

    # Clip data to in/out water times or via good_ens
    ds = utils.clip_ds(ds)

    # rename variables
    ds = ds_rename_vars(ds)

    ds = magvar_correct(ds)

    # should function this
    for var in ds.data_vars:
        ds = qaqc.trim_min(ds, var)
        ds = qaqc.trim_max(ds, var)
        ds = qaqc.trim_min_diff(ds, var)
        ds = qaqc.trim_max_diff(ds, var)
        ds = qaqc.trim_med_diff(ds, var)
        ds = qaqc.trim_med_diff_pct(ds, var)
        ds = qaqc.trim_bad_ens(ds, var)
        ds = qaqc.trim_maxabs_diff_2d(ds, var)
        ds = qaqc.trim_fliers(ds, var)

    # after check for masking vars by other vars
    for var in ds.data_vars:
        ds = qaqc.trim_mask(ds, var)

    # check for drop_vars is config yaml
    if "drop_vars" in ds.attrs:
        ds = qaqc.drop_vars(ds)

    ds = utils.create_z(ds)  # added 7/31/2023

    ds = ds_add_attrs(ds)

    ds = utils.add_standard_names(ds)

    # assign min/max:
    ds = utils.add_min_max(ds)

    ds = utils.add_start_stop_time(ds)

    ds = utils.add_delta_t(ds)

    # add lat/lon coordinates
    ds = utils.ds_add_lat_lon(ds)

    if "vert_dim" in ds.attrs:
        vdim = ds.attrs["vert_dim"]
        attrsbak = ds[vdim].attrs
        # axis attr set for z in utils.create_z so need to del if other than z
        if ds.attrs["vert_dim"] != "z":
            ds[vdim].attrs["axis"] = "Z"
            del ds["z"].attrs["axis"]

    # Write to .nc file
    print("Writing cleaned/trimmed data to .nc file")
    nc_filename = ds.attrs["filename"] + "-a.nc"

    ds.to_netcdf(
        nc_filename, unlimited_dims=["time"], encoding={"time": {"dtype": "i4"}}
    )
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])
    print("Done writing netCDF file", nc_filename)
