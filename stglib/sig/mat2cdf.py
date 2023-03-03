import datetime as dt
import glob
import re
import warnings

import numpy as np
import pandas as pd
import scipy.io
import xarray as xr

from ..aqd import aqdutils
from ..core import utils


def matlab2datetime(matlab_datenum):
    day = dt.datetime.fromordinal(int(matlab_datenum))
    dayfrac = dt.timedelta(days=matlab_datenum % 1) - dt.timedelta(days=366)
    return day + dayfrac


def load_mat_file(filnam):
    mat = utils.loadmat(filnam)
    bindist = (
        mat["Config"]["Burst_BlankingDistance"]
        + mat["Config"]["Burst_CellSize"] / 2
        + mat["Config"]["Burst_CellSize"] * np.arange(mat["Config"]["Burst_NCells"])
    )

    # BurstRawAltimeter
    dsbra = xr.Dataset()
    dsbra["time"] = xr.DataArray(
        [matlab2datetime(x) for x in mat["Data"]["BurstRawAltimeter_Time"]],
        dims="time",
    )
    dsbra["time"] = pd.DatetimeIndex(dsbra["time"])
    dsbra["time"] = pd.DatetimeIndex(dsbra["time"])
    dsbra.attrs["data_type"] = "BurstRawAltimeter"

    # IBurst
    dsi = xr.Dataset()
    dsi["time"] = pd.DatetimeIndex(
        xr.DataArray(
            [matlab2datetime(x) for x in mat["Data"]["IBurst_Time"]],
            dims="time",
        )
    )
    dsi["time"] = pd.DatetimeIndex(dsi["time"])
    dsi["bindist"] = xr.DataArray(bindist, dims="bindist")
    dsi.attrs["data_type"] = "IBurst"

    # Burst
    dsb = xr.Dataset()
    dsb["time"] = pd.DatetimeIndex(
        xr.DataArray(
            [matlab2datetime(x) for x in mat["Data"]["Burst_Time"]],
            dims="time",
        )
    )
    dsb["time"] = pd.DatetimeIndex(dsb["time"])
    dsb["bindist"] = xr.DataArray(bindist, dims="bindist")
    dsb.attrs["data_type"] = "Burst"

    for k in mat["Data"]:
        if "BurstRawAltimeter" in k:
            if "_Time" not in k:
                if mat["Data"][k].ndim == 1:
                    dsbra[k.split("_")[1]] = xr.DataArray(mat["Data"][k], dims="time")
                elif mat["Data"][k].ndim == 2:
                    # print("still need to process", k, mat["Data"][k].shape)
                    pass
        elif "IBurst" in k:
            if "_Time" not in k:
                if mat["Data"][k].ndim == 1:
                    dsi[k.split("_")[1]] = xr.DataArray(mat["Data"][k], dims="time")
                elif mat["Data"][k].ndim == 2:
                    # only checks to see if cells match on first sample
                    if mat["Data"][k].shape[1] == mat["Data"]["IBurst_NCells"][0]:
                        dsi[k.split("_")[1]] = xr.DataArray(
                            mat["Data"][k], dims=["time", "bindist"]
                        )
                else:
                    print("still need to process", k, mat["Data"][k].shape)
        elif re.match("^Burst_", k):
            if "_Time" not in k:
                if mat["Data"][k].ndim == 1:
                    dsb[k.split("_")[1]] = xr.DataArray(mat["Data"][k], dims="time")
                elif mat["Data"][k].ndim == 2:
                    # only checks to see if cells match on first sample
                    if mat["Data"][k].shape[1] == mat["Data"]["Burst_NCells"][0]:
                        dsb[k.split("_")[1]] = xr.DataArray(
                            mat["Data"][k], dims=["time", "bindist"]
                        )
                else:
                    print("still need to process", k, mat["Data"][k].shape)
        else:
            print("missing variable:", k)

    for ds in [dsbra, dsi, dsb]:
        read_config_mat(mat, ds)

    for ds in [dsbra, dsi, dsb]:
        print(ds.attrs["data_type"])
        add_descriptions(mat, ds)

    for ds in [dsbra, dsi, dsb]:
        add_units(mat, ds)

    return dsbra, dsi, dsb


def read_config_mat(mat, ds):
    for k in mat["Config"]:
        if k == "Burst_Beam2xyz":
            ds.attrs[f"SIG{k}"] = str(mat["Config"][k])
        else:
            ds.attrs[f"SIG{k}"] = mat["Config"][k]


def add_descriptions(mat, ds):
    for k in mat["Descriptions"]:
        var = k.split("_")[1]
        if var in ds:
            if "long_name" not in ds[var].attrs:
                ds[var].attrs["long_name"] = mat["Descriptions"][k]
            else:
                warnings.warn(
                    f"long_name already exists for {ds.attrs['data_type']} {var} {k}"
                )


def add_units(mat, ds):
    for k in mat["Units"]:
        var = k.split("_")[1]
        if var in ds:
            if "units" not in ds[var].attrs:
                ds[var].attrs["units"] = mat["Units"][k]
            else:
                warnings.warn(
                    f"units already exists for {ds.attrs['data_type']} {var} {k}"
                )


def mat_to_cdf(metadata):

    basefile = metadata["basefile"]

    if "prefix" in metadata:
        prefix = metadata["prefix"]
    else:
        prefix = ""
    basefile = prefix + basefile

    utils.check_valid_globalatts_metadata(metadata)
    aqdutils.check_valid_config_metadata(metadata)

    dsbras = []
    dsis = []
    dsbs = []

    fildir = "/Users/dnowacki/OneDrive - DOI/stglib/sig1k/"

    for f in glob.glob(f"{basefile}_*.mat"):
        a, b, c = load_mat_file(f)
        filstub = f.split("/")[-1].split(".mat")[0]
        dsbras.append(a)
        dsis.append(b)
        dsbs.append(c)

    dsbra = xr.merge(dsbras)
    dsi = xr.merge(dsis)
    dsb = xr.merge(dsbs)

    # Won't deal with BurstRawAltimeter data for now
    for ds in [dsi, dsb]:
        ds = utils.write_metadata(ds, metadata)
        ds = utils.ensure_cf(ds)

    # Won't deal with BurstRawAltimeter data for now
    # cdf_filename = prefix + dsbra.attrs["filename"] + "bra-raw.cdf"
    # dsbra.to_netcdf(cdf_filename, unlimited_dims=["time"])
    # print(f"Finished writing data to {cdf_filename}")

    cdf_filename = prefix + dsi.attrs["filename"] + "iburst-raw.cdf"
    print(dsi.data_vars)
    dsi.to_netcdf(cdf_filename, unlimited_dims=["time"])
    print(f"Finished writing data to {cdf_filename}")

    cdf_filename = prefix + dsb.attrs["filename"] + "burst-raw.cdf"
    print(dsb.data_vars)
    dsb.to_netcdf(cdf_filename, unlimited_dims=["time"])
    print(f"Finished writing data to {cdf_filename}")
