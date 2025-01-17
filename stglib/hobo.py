import csv
import string

import numpy as np
import pandas as pd
import xarray as xr

from .core import qaqc, utils


def read_hobo(
    filnam, skiprows=1, skipfooter=0, names=["#", "DateTime", "AbsPres_kPa", "Temp_C"]
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

    hobo.set_index("time", inplace=True)

    return xr.Dataset(hobo)


def csv_to_cdf(metadata):
    """
    Process HOBO .csv file to a raw .cdf file
    """

    basefile = metadata["basefile"]

    kwargs = {"skiprows": metadata["skiprows"], "skipfooter": metadata["skipfooter"]}
    if "names" in metadata:
        kwargs["names"] = metadata["names"]
    else:
        kwargs["names"] = get_col_names(basefile + ".csv", metadata)

    ds = read_hobo(basefile + ".csv", **kwargs)

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
    todrop = ["#", "DateTime"]
    return ds.drop([x for x in todrop if x in ds])


def ds_rename_vars(ds):
    # convert some units

    # convert from µS/cm to S/m
    for v in [
        "Conductance_uSpercm",
        "HighRange_Spercm",
        "LowRange_Spercm",
        "SpecificConductance_uSpercm",
    ]:
        if v in ds:
            ds[v].values = ds[v].values / 10000

    # convert from kPa to millibars
    if "AbsPresBarom_kPa" in ds:
        ds["AbsPresBarom_kPa"].values = ds["AbsPresBarom_kPa"].values * 10
        ds = ds.rename({"AbsPresBarom_kPa": "AbsPresBarom_mbar"})

    # convert from kPa to decibars
    if "AbsPres_kPa" in ds:
        ds["AbsPres_kPa"].values = ds["AbsPres_kPa"].values / 10
        ds = ds.rename({"AbsPres_kPa": "AbsPres_dbar"})

    # check to see if logger was deployed as barometer
    if ds.attrs["instrument_type"] == "hwlb":
        # convert from dbar to millibars
        if "AbsPres_dbar" in ds:
            ds["AbsPres_dbar"].values = ds["AbsPres_dbar"].values * 100
            ds = ds.rename({"AbsPres_dbar": "AbsPresBarom_mbar"})
        if "Temp_C" in ds:
            ds = ds.rename({"Temp_C": "Atemp_C"})

    # set up dict of instrument -> EPIC variable names
    varnames = {
        "AbsPres_dbar": "P_1",
        "Temp_C": "T_28",
        "AbsPresBarom_mbar": "BPR_915",
        "SensorDepth_meters": "D_3",
        "Conductance_uSpercm": "C_51",
        "SpecificConductance_uSpercm": "SpC_48",
        "Salinity_ppt": "S_41",
        "DOPercentSat_percent": "OST_62",
        "DOconc_mgperL": "DO",
        "DOAdjConc_mgperL": "DO_Adj",
        "Atemp_C": "T_21",
        "HighRange_Spercm": "C_51_hi",
        "LowRange_Spercm": "C_51_lo",
    }

    # check to make sure they exist before trying to rename
    newvars = {}
    for k in varnames:
        if k in ds:
            newvars[k] = varnames[k]

    # drop unneeded vars
    todrop = ["FullRange_uSpercm"]
    ds = ds.drop([x for x in todrop if x in ds])

    return ds.rename(newvars)


def ds_add_attrs(ds):
    # Update attributes for EPIC and STG compliance
    ds = utils.ds_coord_no_fillvalue(ds)

    ds["time"].attrs.update(
        {"standard_name": "time", "axis": "T", "long_name": "time (UTC)"}
    )

    # Some legacy code leave for now -part for conductivity
    if "condlo_uScm" in ds:
        ds = ds.rename({"condlo_uScm": "SpC_48_lo"})

        ds["SpC_48_lo"].attrs.update(
            {
                "units": "uS/cm",
                "long_name": "Conductivity",
                "comment": "Temperature compensated to 25 °C; low range",
                "epic_code": 48,
                "standard_name": "sea_water_electrical_conductivity_at_reference_temperature",
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
                "standard_name": "sea_water_electrical_conductivity_at_reference_temperature",
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
    # End of legacy code section

    if "T_28" in ds:
        ds["T_28"].attrs.update(
            {
                "units": "degree_C",
                "long_name": "Temperature",
                "epic_code": 28,
                "standard_name": "sea_water_temperature",
            }
        )

    for v in ["C_51", "C_51_hi", "C_51_lo"]:
        if v in ds:
            long_name = "Conductivity"
            if v == "C_51_hi":
                long_name = long_name + ", high range (0.5 to 5.5 S/m)"
            elif v == "C_51_lo":
                long_name = long_name + ", low range (0.01 to 1 S/m)"
            ds[v].attrs.update(
                {
                    "units": "S/m",
                    "long_name": long_name,
                    "epic_code": 51,
                    "standard_name": "sea_water_electrical_conductivity",
                }
            )

    if "SpC_48" in ds:
        ds["SpC_48"].attrs.update(
            {
                "units": "S/m",
                "long_name": "Specific Conductivity",
                "comment": "Temperature compensated to 25 C",
                "epic_code": 48,
                "standard_name": "sea_water_electrical_conductivity",
            }
        )

    if "S_41" in ds:
        ds["S_41"].attrs.update(
            {
                "units": "1",
                "long_name": "Salinity, PSU",
                "comments": "Practical salinity units (PSU)",
                "epic_code": 41,
                "standard_name": "sea_water_practical_salinity",
            }
        )

    if "OST_62" in ds:
        ds["OST_62"].attrs.update(
            {
                "units": "percent",
                "long_name": "Oxygen percent saturation",
                "epic_code": 62,
                "standard_name": "fractional_saturation_of_oxygen_in_sea_water",
            }
        )

    if "DO" in ds:
        ds["DO"].attrs.update(
            {
                "units": "mg/L",
                "long_name": "Dissolved oxygen",
                "standard_name": "mass_concentration_of_oxygen_in_sea_water",
            }
        )

    if "DO_Adj" in ds:
        ds["DO"].values = ds["DO_Adj"].values
        if "DO_note" in ds.attrs:
            # ds = utils.insert_note(ds, "DO", ds.attrs["DO_note"] + " ")
            ds["DO"].attrs.update({"note": ds.attrs["DO_note"]})
        else:
            ds["DO"].attrs.update(
                {
                    "note": "Using adjusted DO concentration",
                }
            )

        ds = ds.drop_vars("DO_Adj")

    if "P_1" in ds:
        ds["P_1"].attrs.update(
            {
                "units": "dbar",
                "long_name": "Uncorrected pressure",
                "epic_code": 1,
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
            # ds = utils.insert_note(ds, "P_1ac", ds.attrs["P_1ac_note"] + " ")
            ds["P_1ac"].attrs.update({"note": ds.attrs["P_1ac"]})

    if "D_3" in ds:
        ds["D_3"].attrs.update(
            {
                "units": f"{ds.depth.attrs['units']}",
                "long_name": "depth below sea surface",
                "standard_name": "depth",
                "positive": f"{ds.depth.attrs['positive']}",
            }
        )
        if "D_3_note" in ds.attrs:
            # ds = utils.insert_note(ds, "D_3", ds.attrs["D_3_note"] + " ")
            ds["D_3"].attrs.update({"note": ds.attrs["D_3_note"]})

    if "BPR_915" in ds:
        ds["BPR_915"].attrs.update(
            {
                "units": "mbar",
                "long_name": "Barometric pressure",
                "epic_code": 915,
                "standard_name": "air_pressure",
            }
        )

    if "T_21" in ds:
        ds["T_21"].attrs.update(
            {
                "units": "degree_C",
                "long_name": "Air Temperature",
                "epic_code": 21,
                "standard_name": "air_temperature",
            }
        )

    return ds


def get_serial_number(filnam):
    """get the serial number of the instrument"""
    with open(filnam) as f:
        f.readline()
        line2 = f.readline()
        sn = line2.split("LGR S/N: ")[1].split(",")[0]

        return sn


def strip_non_printable(strin):
    """Returns the string without non printable characters"""
    printable = set(string.printable)
    return "".join(filter(lambda x: x in printable, strin))


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
            dcols[k] = dcols[k].replace("/", "per")

        # then strip non-ascii characters
        dcols[k] = strip_non_printable(dcols[k])

    names = []
    for k in dcols:
        if k == "#" or k == "DateTime":
            names.append(k)
        else:
            names.append(k + "_" + dcols[k])

    return names


def cdf_to_nc(cdf_filename, atmpres=False):
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

    if atmpres:
        ds = utils.atmos_correct(ds, atmpres)

    # QAQC
    ds = qaqc.call_qaqc(ds)

    ds = utils.create_z(ds)  # added 7/31/2023
    ds = utils.create_water_level_var(ds)
    ds = utils.create_filtered_water_level_var(ds)
    ds = ds_add_attrs(ds)
    ds = utils.add_min_max(ds)
    ds = utils.add_start_stop_time(ds)
    ds = utils.add_delta_t(ds)
    ds = utils.ds_add_lat_lon(ds)

    if "vert_dim" in ds.attrs:
        vdim = ds.attrs["vert_dim"]
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
