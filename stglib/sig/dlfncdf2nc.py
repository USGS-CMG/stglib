import pkgutil
import time

import xarray as xr

import dolfyn
import xmltodict

from ..aqd import aqdutils
from ..core import utils


def cdf_to_nc(cdf_filename, atmpres=False):
    """
    Load a "raw" .cdf file and generate a processed .nc file
    """
    # TODO: Add atmospheric pressure offset
    print(f"Loading {cdf_filename}")
    start_time = time.time()
    ds = dolfyn.load(cdf_filename)
    end_time = time.time()
    print(f"Finished loading {cdf_filename} in {end_time-start_time:.1f} seconds")

    ds = utils.create_nominal_instrument_depth(ds)

    # Transform ("rotate") into ENU coordinates
    dolfyn.rotate2(ds, "earth")

    # Clip data to in/out water times or via good_ens
    # Need to clip data after coord transform when using dolfyn
    ds = utils.clip_ds(ds)

    # Create separate vel variables first
    ds["U"] = ds["vel"].sel(dir="E")
    ds["V"] = ds["vel"].sel(dir="N")
    ds["W1"] = ds["vel"].sel(dir="U1")
    ds["W2"] = ds["vel"].sel(dir="U2")

    ds = aqdutils.magvar_correct(ds)

    # Rename DataArrays for EPIC compliance
    ds = aqdutils.ds_rename(ds)

    # Add EPIC and CMG attributes
    # ds = aqdutils.ds_add_attrs(ds, inst_type="SIG")

    # Add start_time and stop_time attrs
    ds = utils.add_start_stop_time(ds)

    # Add history showing file used
    ds = utils.add_history(ds)

    ds = clean_dolfyn_standard_names(ds)

    ds = utils.add_standard_names(ds)

    ds = fix_dolfyn_encoding(ds)

    # split up into multiple files and write them separately
    ds_b5 = ds.copy()
    ds_echo = ds.copy()

    todrop = []
    for v in ds.data_vars:
        if "_b5" in v or "_echo" in v:
            todrop.append(v)

    ds = ds.drop_vars(todrop)
    ds = drop_unused_dims(ds)

    todrop = []
    for v in ds_b5.data_vars:
        if "_b5" not in v:
            todrop.append(v)

    ds_b5 = ds_b5.drop_vars(todrop)
    ds_b5 = drop_unused_dims(ds_b5)

    todrop = []
    for v in ds_echo.data_vars:
        if "_echo" not in v:
            todrop.append(v)

    ds_echo = ds_echo.drop_vars(todrop)
    ds_echo = drop_unused_dims(ds_echo)

    if "prefix" in ds.attrs:
        nc_filename = ds.attrs["prefix"] + ds.attrs["filename"] + "-a.nc"
    else:
        nc_filename = ds.attrs["filename"] + "-a.nc"

    for datatype, dsout in zip(["", "_b5", "_echo"], [ds, ds_b5, ds_echo]):
        if datatype == "":
            nc_out = nc_filename
        elif datatype == "_b5":
            nc_out = nc_filename[:-5] + "_b5-a.nc"
        elif datatype == "_echo":
            nc_out = nc_filename[:-5] + "_echo-a.nc"

        dolfyn.save(dsout, nc_out)
        utils.check_compliance(nc_out, conventions=ds.attrs["Conventions"])

        print("Done writing netCDF file", nc_out)

    return ds


def drop_unused_dims(ds):
    """only keep dims that will be in the final files"""
    thedims = []
    for v in ds.data_vars:
        for x in ds[v].dims:
            thedims.append(x)

    for x in ds.dims:
        if x not in thedims:
            ds = ds.drop_vars(x)

    return ds


def clean_dolfyn_standard_names(ds):
    """remove non-compliant standard_names set by dolfyn"""
    data = pkgutil.get_data(__name__, "../data/cf-standard-name-table.xml")
    doc = xmltodict.parse(data)

    entries = doc["standard_name_table"]["entry"]
    allnames = [x["@id"] for x in entries]

    for v in ds.data_vars:
        if "standard_name" in ds[v].attrs:
            if ds[v].attrs["standard_name"] not in allnames:
                del ds[v].attrs["standard_name"]

    return ds


def fix_dolfyn_encoding(ds):
    """ensure we don't set dtypes of int64 for CF compliance"""
    for var in ds.dims:
        if ds[var].dtype == "int64":
            if ds[var].max() > 2**31 - 1 or ds[var].min() < -(2**31):
                print(
                    f"warning {var} may be too big to fit in int32: min {ds[var].min().values}, max {ds[var].max().values}"
                )
            ds[var].encoding["dtype"] = "int32"

    return ds
