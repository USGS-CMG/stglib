import pandas as pd
import xarray as xr

from .core import qaqc, utils


def read_asc(filnam, skiprows=51, encoding="utf-8"):
    """Read data from an SBE 37 MicroCAT .asc file into an xarray
    Dataset.

    Parameters
    ----------
    filnam : string
        The filename
    skiprows : int, optional
        How many header rows to skip. Default 51
    encoding : string, optional
        File encoding. Default 'utf-8'
    Returns
    -------
    xarray.Dataset
        An xarray Dataset of the MicroCAT data
    """

    df = pd.read_csv(
        filnam,
        skiprows=skiprows,
        header=None,
        names=["Temp", "Cond", "Sal", "Date", "Time"],
        parse_dates={"time": ["Date", "Time"]},
        encoding=encoding,
        index_col=False,
    )
    df.set_index("time", inplace=True)
    mc = df.to_xarray()
    return mc


def read_asc_header(filnam):
    """Read header from an SBE 37 MicroCAT .asc file and add to metadata"""

    header = {}

    with open(filnam) as file:
        for line in file:
            col = line.split()
            if "*" not in line:
                break
            elif "SERIAL NO." in line:
                header["serial_number"] = col[6]
                header["instrument_type"] = col[1]
            elif "sample interval" in line:
                header["sample_interval"] = col[4]

    return header


def asc_to_cdf(metadata):
    """
    Load a raw .asc file and generate a .cdf file
    """
    basefile = metadata["basefile"]

    ds = read_asc(basefile + ".asc", skiprows=metadata["skiprows"])

    header = read_asc_header(basefile + ".asc")

    metadata.update(header)

    metadata.pop("skiprows")

    ds = utils.write_metadata(ds, metadata)

    ds = utils.ensure_cf(ds)

    # configure file
    cdf_filename = ds.attrs["filename"] + "-raw.cdf"

    ds.to_netcdf(cdf_filename, unlimited_dims=["time"])

    print(f"Finished writing data to {cdf_filename}")

    return ds


def cdf_to_nc(cdf_filename):
    """
    Load a raw .cdf file and generate a processed .nc file
    """

    # Load raw .cdf data
    ds = xr.open_dataset(cdf_filename)

    # remove units in case we change and we can use larger time steps
    ds.time.encoding.pop("units")

    # Rename variables to CF compliant names
    ds = ds_rename_vars(ds)

    # Add attributes
    ds = ds_add_attrs(ds)

    # QAQC
    ds = qaqc.call_qaqc(ds)

    # Run utilities
    ds = utils.create_z(ds)
    ds = utils.clip_ds(ds)
    ds = utils.ds_add_lat_lon(ds)
    ds = utils.create_nominal_instrument_depth(ds)
    ds = utils.add_start_stop_time(ds)
    ds = utils.add_min_max(ds)
    ds = utils.add_delta_t(ds)
    ds = utils.ds_coord_no_fillvalue(ds)

    # Write to .nc file
    print("Writing cleaned/trimmed data to .nc file")
    nc_filename = ds.attrs["filename"] + "-a.nc"

    ds.to_netcdf(nc_filename, unlimited_dims=["time"])
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])

    print(f"Done writing netCDF file {nc_filename}")


def ds_rename_vars(ds):
    """
    Rename variables to be CF compliant
    """
    varnames = {"Temp": "T_28", "Cond": "C_51", "Sal": "S_41"}

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

    if "C_51" in ds:
        ds["C_51"].attrs.update(
            {
                "units": "S/m",
                "long_name": "Conductivity",
                "standard_name": "sea_water_electrical_conductivity",
            }
        )

    if "S_41" in ds:
        ds["S_41"].attrs.update(
            {
                "units": "1",
                "long_name": "Salinity, PSU",
                "comments": "Practical salinity units (PSU)",
                "standard_name": "sea_water_practical_salinity",
            }
        )

    return ds
