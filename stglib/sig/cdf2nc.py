import time

import dolfyn

from ..aqd import aqdutils
from ..core import utils


def cdf_to_nc(cdf_filename, atmpres=False):
    """
    Load a "raw" .cdf file and generate a processed .nc file
    """

    # TODO: Add atmospheric pressure offset
    print(f"Loading {cdf_filename[0]}")
    start_time = time.time()
    ds = dolfyn.load(cdf_filename[0])
    end_time = time.time()
    print(f"Finished loading {cdf_filename[0]} in {end_time-start_time:.1f} seconds")

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

    ds = set_dolfyn_data_types(ds)

    # ds = utils.add_standard_names(ds)

    # split up into multiple files:
    ds_b5 = ds.copy()
    ds_echo = ds.copy()

    todrop = []
    for v in ds.data_vars:
        if "_b5" in v:
            todrop.append(v)

    ds = ds.drop_vars(todrop)

    todrop = []
    for v in ds_b5.data_vars:
        if "_b5" not in v:
            todrop.append(v)

    ds_b5 = ds_b5.drop_vars(todrop)

    todrop = []
    for v in ds_echo.data_vars:
        if "_echo" not in v:
            todrop.append(v)

    ds_echo = ds_echo.drop_vars(todrop)

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

        if datatype != "":
            dolfyn.save(dsout, nc_out, compression=True)
            utils.check_compliance(nc_out, conventions=ds.attrs["Conventions"])

            print("Done writing netCDF file", nc_filename)

    return ds


def clean_dolfyn_standard_names(ds):
    import xmltodict

    with open("cf-standard-name-table.xml", encoding="utf-8") as fd:
        doc = xmltodict.parse(fd.read())

    entries = doc["standard_name_table"]["entry"]
    allnames = [x["@id"] for x in entries]

    for v in ds.data_vars:
        if "standard_name" in ds[v].attrs:
            if ds[v].attrs["standard_name"] not in allnames:
                print("removing", v, "standard_name")
                del ds[v].attrs["standard_name"]
                print(v, f"{ds[v].attrs}")

    return ds


def set_dolfyn_data_types(ds):
    """make datatypes for time, etc not be int64"""
    for d in ds.dims:
        print(d)
        print(ds[d].dtype)
        print(ds[d].max())
    return ds