from __future__ import division, print_function

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
    hobo["time"] = pd.to_datetime(hobo["datetime"])
    if "abspres_kPa" in hobo:
        hobo["abspres_dbar"] = hobo["abspres_kPa"] / 10
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
    try:
        ds = read_hobo(basefile + ".csv", **kwargs)
    except UnicodeDecodeError:
        # try reading as Mac OS Western for old versions of Mac Excel
        ds = read_hobo(basefile + ".csv", encoding="mac-roman", **kwargs)

    metadata.pop("skiprows")
    metadata.pop("skipfooter")

    # write out metadata first, then deal exclusively with xarray attrs
    ds = utils.write_metadata(ds, metadata)

    del metadata

    ds = utils.shift_time(ds, 0)

    if not utils.is_cf(ds):
        ds = utils.create_epic_times(ds)

    ds = drop_vars(ds)

    ds.attrs["serial_number"] = get_serial_number(basefile + ".csv")

    # configure file
    cdf_filename = ds.attrs["filename"] + "-raw.cdf"

    ds.to_netcdf(cdf_filename, unlimited_dims=["time"])

    print("Finished writing data to %s" % cdf_filename)

    return ds


def drop_vars(ds):
    todrop = ["#", "datetime", "abspres_kPa"]
    return ds.drop([x for x in todrop if x in ds])


def ds_add_attrs(ds):

    # Update attributes for EPIC and STG compliance
    ds = utils.ds_coord_no_fillvalue(ds)

    ds["time"].attrs.update({"standard_name": "time", "axis": "T"})

    if not utils.is_cf(ds):
        ds["epic_time"].attrs.update(
            {"units": "True Julian Day", "type": "EVEN", "epic_code": 624}
        )

        ds["epic_time2"].attrs.update(
            {"units": "msec since 0:00 GMT", "type": "EVEN", "epic_code": 624}
        )

    if "abspres_dbar" in ds:
        ds = ds.rename({"abspres_dbar": "BPR_915"})

        # convert decibar to millibar
        ds["BPR_915"] = ds["BPR_915"] * 100

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

        ds["SpC_48_lo"].values = (
            ds["SpC_48_lo"].values / 1000
        )  # convert from µS/cm to mS/cm

        ds["SpC_48_lo"].attrs.update(
            {
                "units": "mS/cm",
                "long_name": "Conductivity",
                "comment": "Temperature compensated to 25 °C; low range",
                "epic_code": 48,
                "standard_name": "sea_water_electrical_conductivity",
            }
        )

        ds["S_41_lo"] = utils.salinity_from_spcon(ds["SpC_48_lo"])

        ds["S_41_lo"].attrs.update(
            {
                "units": "1e-3",
                "long_name": "Salinity; low range",
                "epic_code": 41,
                "standard_name": "sea_water_salinity",
            }
        )

    if "condhi_uScm" in ds:
        ds = ds.rename({"condhi_uScm": "SpC_48_hi"})

        ds["SpC_48_hi"].values = (
            ds["SpC_48_hi"].values / 1000
        )  # convert from µS/cm to mS/cm

        ds["SpC_48_hi"].attrs.update(
            {
                "units": "mS/cm",
                "long_name": "Conductivity",
                "comment": "Temperature compensated to 25 °C; high range",
                "epic_code": 48,
                "standard_name": "sea_water_electrical_conductivity",
            }
        )

        ds["S_41_hi"] = utils.salinity_from_spcon(ds["SpC_48_hi"])

        ds["S_41_hi"].attrs.update(
            {
                "units": "1e-3",
                "long_name": "Salinity; high range",
                "epic_code": 41,
                "standard_name": "sea_water_salinity",
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
        var.encoding["_FillValue"] = 1e35

    for var in ds.variables:
        if (var not in ds.coords) and ("time" not in var):
            add_attributes(ds[var], ds.attrs)

    ds.attrs["COMPOSITE"] = np.int32(0)

    return ds


def get_serial_number(filnam):
    """get the serial number of the instrument"""

    with open(filnam) as f:
        f.readline()
        line2 = f.readline()
        sn = line2.find("LGR S/N: ")
        # these are the indices of the serial number
        return line2[sn + 9 : sn + 17]


def cdf_to_nc(cdf_filename):
    """
    Load a "raw" .cdf file and generate a processed .nc file
    """

    # Load raw .cdf data
    ds = xr.open_dataset(cdf_filename)
    ds.time.encoding.pop(
        "units"
    )  # remove units in case we change and we can use larger time steps

    # Clip data to in/out water times or via good_ens
    ds = utils.clip_ds(ds)

    # assign min/max:
    ds = utils.add_min_max(ds)

    ds = utils.add_start_stop_time(ds)

    if not utils.is_cf(ds):
        ds = utils.create_epic_times(ds)

    ds = utils.add_delta_t(ds)

    # add lat/lon coordinates
    ds = utils.ds_add_lat_lon(ds)

    ds = ds_add_attrs(ds)

    ds = utils.no_p_create_depth(ds)

    # add lat/lon coordinates to each variable
    for var in ds.variables:
        if (var not in ds.coords) and ("time" not in var):
            ds = utils.add_lat_lon(ds, var)
            ds = utils.no_p_add_depth(ds, var)
            # cast as float32
            ds = utils.set_var_dtype(ds, var)

    ds = utils.rename_time(ds)

    # Write to .nc file
    print("Writing cleaned/trimmed data to .nc file")
    nc_filename = ds.attrs["filename"] + "-a.nc"

    ds.to_netcdf(nc_filename, unlimited_dims=["time"])
    print("Done writing netCDF file", nc_filename)
