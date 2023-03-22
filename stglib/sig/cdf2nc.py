import math

import dolfyn
import numpy as np
import xarray as xr
from tqdm import tqdm

from ..aqd import aqdutils
from ..core import qaqc, utils


def cdf_to_nc(cdf_filename, atmpres=False):
    """
    Load a "raw" .cdf file and generate a processed .nc file
    """

    # TODO: Add atmospheric pressure offset
    ds = dolfyn.load(cdf_filename[0])

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

    # ds = utils.add_standard_names(ds)

    if "prefix" in ds.attrs:
        nc_filename = ds.attrs["prefix"] + ds.attrs["filename"] + "-a.nc"
    else:
        nc_filename = ds.attrs["filename"] + "-a.nc"

    dolfyn.save(ds, nc_filename, compression=True)
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])

    print("Done writing netCDF file", nc_filename)

    return ds
