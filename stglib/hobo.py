import csv

import numpy as np
import pandas as pd
import xarray as xr

from .core import utils


def read_hobo(
    filnam, skiprows=1, skipfooter=0, names=["#", "datetime", "abspres_kPa", "temp_C"]
):
    """Read data from an Onset HOBO pressure sensor .csv file into an xarray
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
        An xarray Dataset of the HOBO data
    """
    hobo = pd.read_csv(
        filnam,
        usecols=np.arange(len(names)),
        names=names,
        engine="python",
        skiprows=skiprows,
        skipfooter=skipfooter,
    )

    hobo["time"] = pd.to_datetime(hobo["DateTime"])
    if "AbsPres_kPa" in hobo:
        hobo["AbsPres_dbar"] = hobo["AbsPres_kPa"] / 10
    hobo.set_index("time", inplace=True)

    return xr.Dataset(hobo)


def csv_to_cdf(metadata):
    """
    Process HOBO .csv file to a raw .cdf file
    """

    basefile = metadata["basefile"]

    names = get_col_names(basefile + ".csv", metadata)

    kwargs = {"skiprows": metadata["skiprows"], "skipfooter": metadata["skipfooter"]}
    # if "names" in metadata:
    #    kwargs["names"] = metadata["names"]
    kwargs["names"] = names
    try:
        ds = read_hobo(basefile + ".csv", **kwargs)
    except UnicodeDecodeError:
        # try reading as Mac OS Western for old versions of Mac Excel
        ds = read_hobo(basefile + ".csv", encoding="mac_roman", **kwargs)

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

    ds.attrs["serial_number"] = get_serial_number(basefile + ".csv")

    # configure file
    cdf_filename = ds.attrs["filename"] + "-raw.cdf"

    ds.to_netcdf(cdf_filename, unlimited_dims=["time"])

    print("Finished writing data to %s" % cdf_filename)

    return ds


def drop_vars(ds):
    todrop = ["#", "DateTime", "AbsPres_kPa"]
    return ds.drop([x for x in todrop if x in ds])


def ds_add_attrs(ds):
    # Update attributes for EPIC and STG compliance
    ds = utils.ds_coord_no_fillvalue(ds)

    ds["time"].attrs.update(
        {"standard_name": "time", "axis": "T", "long_name": "time (UTC)"}
    )

    if "abspres_dbar" in ds:
        ds = ds.rename({"abspres_dbar": "P_1"})

        ds["BPR_915"].attrs.update(
            {"units": "mbar", "long_name": "Barometric pressure", "epic_code": 915}
        )

    if "baropres_kPa" in ds:
        ds = ds.rename({"baropres_kPa": "BPR_915"})

        # convert kPa to millibar
        ds["BPR_915"] = ds["BPR_915"] * 10

        ds["BPR_915"].attrs.update(
            {"units": "mbar", "long_name": "Barometric pressure", "epic_code": 915}
        )

    if "temp_C" in ds:
        ds = ds.rename({"temp_C": "T_28"})

        ds["T_28"].attrs.update(
            {
                "units": "C",
                "long_name": "Temperature",
                "epic_code": 28,
                "standard_name": "sea_water_temperature",
            }
        )

    if "condlo_uScm" in ds:
        ds = ds.rename({"condlo_uScm": "SpC_48_lo"})

        ds["SpC_48_lo"].attrs.update(
            {
                "units": "uS/cm",
                "long_name": "Conductivity",
                "comment": "Temperature compensated to 25 °C; low range",
                "epic_code": 48,
                "standard_name": "sea_water_electrical_conductivity",
            }
        )

        ds["S_41_lo"] = utils.salinity_from_spcon(ds["SpC_48_lo"])

        ds["S_41_lo"].attrs.update(
            {
                "units": "1",
                "long_name": "Salinity; low range, PSU",
                "epic_code": 41,
                "standard_name": "sea_water_practical_salinity",
            }
        )

    if "condhi_uScm" in ds:
        ds = ds.rename({"condhi_uScm": "SpC_48_hi"})

        ds["SpC_48_hi"].attrs.update(
            {
                "units": "uS/cm",
                "long_name": "Conductivity",
                "comment": "Temperature compensated to 25 °C; high range",
                "epic_code": 48,
                "standard_name": "sea_water_electrical_conductivity",
            }
        )

        ds["S_41_hi"] = utils.salinity_from_spcon(ds["SpC_48_hi"])

        ds["S_41_hi"].attrs.update(
            {
                "units": "1",
                "long_name": "Salinity; high range, PSU",
                "epic_code": 41,
                "standard_name": "sea_water_practical_salinity",
            }
        )

    def add_attributes(var, dsattrs):
        var.attrs.update(
            {
                "initial_instrument_height": dsattrs["initial_instrument_height"],
                # 'nominal_instrument_depth': dsattrs['nominal_instrument_depth'],
                "height_depth_units": "m",
            }
        )

    for var in ds.variables:
        if (var not in ds.coords) and ("time" not in var):
            add_attributes(ds[var], ds.attrs)

    return ds


def get_serial_number(filnam):
    """get the serial number of the instrument"""
    with open(filnam) as f:
        f.readline()
        line2 = f.readline()
        sn = line2.split("LGR S/N: ")[1].split(",")[0]

        return sn


def get_col_names(filnam, metadata):
    """get column names and column units from instrument input data file"""

    with open(filnam) as f:
        rdr = csv.reader(f)
        # check to see if first value is "#"
        hdrline = next(rdr)
        while hdrline[0] != "#":
            hdrline = next(rdr)

    collist = [x.split(" (")[0] for x in hdrline]

    colnames = []
    colunits = []
    for x in collist:
        # spl = x.split(",")
        spl = x.replace(" ", "").split(",")
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
            dcols[k] = dcols[k].replace("/", "per")
        if "uSpercm" in dcols[k]:
            dcols[k] = "uSpercm"

    names = []
    for k in dcols:
        if k == "#" or k == "DateTime":
            names.append(k)
        else:
            names.append(k + "_" + dcols[k])

    return names


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

    # assign min/max:
    ds = utils.add_min_max(ds)

    ds = utils.add_start_stop_time(ds)

    ds = utils.add_delta_t(ds)

    # add lat/lon coordinates
    ds = utils.ds_add_lat_lon(ds)

    ds = ds_add_attrs(ds)

    ds = utils.no_p_create_depth(ds)

    # add lat/lon coordinates to each variable
    for var in ds.variables:
        if (var not in ds.coords) and ("time" not in var):
            # ds = utils.add_lat_lon(ds, var)
            ds = utils.no_p_add_depth(ds, var)
            # cast as float32
            ds = utils.set_var_dtype(ds, var)

    # Write to .nc file
    print("Writing cleaned/trimmed data to .nc file")
    nc_filename = ds.attrs["filename"] + "-a.nc"

    ds.to_netcdf(nc_filename, unlimited_dims=["time"])
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])
    print("Done writing netCDF file", nc_filename)
