import datetime as dt
import glob
import os
import re
import time

import numpy as np
import pandas as pd
import xarray as xr
from dask.diagnostics import ProgressBar
from joblib import Parallel, delayed

from ..aqd import aqdutils
from ..core import utils


def matlab2datetime(matlab_datenum):
    day = dt.datetime.fromordinal(int(matlab_datenum))
    dayfrac = dt.timedelta(days=matlab_datenum % 1) - dt.timedelta(days=366)
    return day + dayfrac


def load_mat_file(filnam):
    mat = utils.loadmat(filnam)

    if (
        mat["Config"]["Burst_HighResolution"] == "False"
        or mat["Config"]["Burst_HighResolution5"] == "False"
    ):
        bindist = (
            mat["Config"]["Burst_BlankingDistance"]
            + mat["Config"]["Burst_CellSize"] / 2
            + mat["Config"]["Burst_CellSize"] * np.arange(mat["Config"]["Burst_NCells"])
        )

    if (
        mat["Config"]["Burst_HighResolution"] == "True"
        or mat["Config"]["Burst_HighResolution5"] == "True"
    ):
        bindistHR = (
            mat["Config"]["BurstHR_BlankingDistance"]
            + mat["Config"]["BurstHR_CellSize"] / 2
            + mat["Config"]["BurstHR_CellSize"]
            * np.arange(mat["Config"]["BurstHR_NCells"])
        )

    if mat["Config"]["Burst_EchoSounder"] == "True":
        bindistECHO = (
            mat["Config"]["EchoSounder_BlankingDistance"]
            + mat["Config"]["EchoSounder_CellSize"] / 2
            + mat["Config"]["EchoSounder_CellSize"]
            * np.arange(mat["Config"]["EchoSounder_NCells"])
        )

    if mat["Config"]["Plan_AverageEnabled"] == "True":
        bindistAVG = (
            mat["Config"]["Average_BlankingDistance"]
            + mat["Config"]["Average_CellSize"] / 2
            + mat["Config"]["Average_CellSize"]
            * np.arange(mat["Config"]["Average_NCells"])
        )

    if "Alt_Plan_AverageEnabled" in mat["Config"]:
        if mat["Config"]["Alt_Plan_AverageEnabled"] == "True":
            bindistAltAVG = (
                mat["Config"]["Alt_Average_BlankingDistance"]
                + mat["Config"]["Alt_Average_CellSize"] / 2
                + mat["Config"]["Alt_Average_CellSize"]
                * np.arange(mat["Config"]["Alt_Average_NCells"])
            )

    ds_dict = {}  # initialize dictionary for xarray datasets
    if (
        mat["Config"]["Plan_BurstEnabled"] == "True"
        and mat["Config"]["Burst_RawAltimeter"] == 1
        and mat["Config"]["Burst_Altimeter"] == "True"
    ):
        # BurstRawAltimeter
        if "BurstRawAltimeter_Time" in mat["Data"]:
            dsbra = xr.Dataset()
            dsbra["time"] = xr.DataArray(
                [matlab2datetime(x) for x in mat["Data"]["BurstRawAltimeter_Time"]],
                dims="time",
            )
            dsbra["time"] = pd.DatetimeIndex(dsbra["time"])
            dsbra["time"] = pd.DatetimeIndex(dsbra["time"])
            dsbra.attrs["data_type"] = "BurstRawAltimeter"
            ds_dict["dsbra"] = dsbra

    if (
        mat["Config"]["Plan_BurstEnabled"] == "True"
        and mat["Config"]["Burst_NBeams"] == 5
    ):
        # IBurst
        if (
            mat["Config"]["Burst_HighResolution5"] == "True"
            and "IBurstHR_Time" in mat["Data"]
        ):
            dsi = xr.Dataset()
            dsi["time"] = pd.DatetimeIndex(
                xr.DataArray(
                    [matlab2datetime(x) for x in mat["Data"]["IBurstHR_Time"]],
                    dims="time",
                )
            )
            dsi["time"] = pd.DatetimeIndex(dsi["time"])
            dsi["bindist"] = xr.DataArray(bindistHR, dims="bindist")
            dsi.attrs["data_type"] = "IBurstHR"
            ds_dict["dsi"] = dsi

        elif "IBurst_Time" in mat["Data"]:
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
            ds_dict["dsi"] = dsi

    if mat["Config"]["Plan_BurstEnabled"] == "True":
        # Burst
        if (
            mat["Config"]["Burst_HighResolution"] == "True"
            and "BurstHR_Time" in mat["Data"]
        ):
            dsb = xr.Dataset()
            dsb["time"] = pd.DatetimeIndex(
                xr.DataArray(
                    [matlab2datetime(x) for x in mat["Data"]["BurstHR_Time"]],
                    dims="time",
                )
            )
            dsb["time"] = pd.DatetimeIndex(dsb["time"])
            dsb["bindist"] = xr.DataArray(bindistHR, dims="bindist")
            dsb.attrs["data_type"] = "BurstHR"
            ds_dict["dsb"] = dsb

        elif "Burst_Time" in mat["Data"]:
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
            ds_dict["dsb"] = dsb

    if (
        mat["Config"]["Plan_BurstEnabled"] == "True"
        and mat["Config"]["Burst_EchoSounder"] == "True"
    ):
        # echo1 data - only handling echo1 data to start
        if "EchoSounder_Frequency1" in mat["Config"]:
            freq1 = mat["Config"]["EchoSounder_Frequency1"]
            if f"Echo1Bin1_{freq1}kHz_Time" in mat["Data"]:
                dse1 = xr.Dataset()
                dse1["time"] = pd.DatetimeIndex(
                    xr.DataArray(
                        [
                            matlab2datetime(x)
                            for x in mat["Data"][f"Echo1Bin1_{freq1}kHz_Time"]
                        ],
                        dims="time",
                    )
                )
                dse1["time"] = pd.DatetimeIndex(dse1["time"])
                dse1["bindist"] = xr.DataArray(bindistECHO, dims="bindist")
                dse1.attrs["data_type"] = "EchoSounder"
                ds_dict["dse1"] = dse1

    if mat["Config"]["Plan_AverageEnabled"] == "True":
        if "Average_Time" in mat["Data"]:
            # Average
            dsa = xr.Dataset()
            dsa["time"] = pd.DatetimeIndex(
                xr.DataArray(
                    [matlab2datetime(x) for x in mat["Data"]["Average_Time"]],
                    dims="time",
                )
            )
            dsa["time"] = pd.DatetimeIndex(dsa["time"])
            dsa["bindist"] = xr.DataArray(bindistAVG, dims="bindist")
            dsa.attrs["data_type"] = "Average"
            ds_dict["dsa"] = dsa

    if "Alt_Plan_AverageEnabled" in mat["Config"]:
        if mat["Config"]["Alt_Plan_AverageEnabled"] == "True":
            if "Alt_Average_Time" in mat["Data"]:
                # Alt Average
                dsalt = xr.Dataset()
                dsalt["time"] = pd.DatetimeIndex(
                    xr.DataArray(
                        [matlab2datetime(x) for x in mat["Data"]["Alt_Average_Time"]],
                        dims="time",
                    )
                )
                dsalt["time"] = pd.DatetimeIndex(dsalt["time"])
                dsalt["bindist"] = xr.DataArray(bindistAltAVG, dims="bindist")
                dsalt.attrs["data_type"] = "Alt_Average"
                ds_dict["dsalt"] = dsalt

    for k in mat["Data"]:
        if "BurstRawAltimeter" in k:
            if "_Time" not in k:
                if mat["Data"][k].ndim == 1:
                    dsbra[k.split("_")[1]] = xr.DataArray(mat["Data"][k], dims="time")
                elif mat["Data"][k].ndim == 2:
                    # print("still need to process", k, mat["Data"][k].shape)
                    pass
        elif "IBurst" in k or "IBurstHR" in k:
            if "_Time" not in k:
                if mat["Data"][k].ndim == 1:
                    dsi[k.split("_")[1]] = xr.DataArray(mat["Data"][k], dims="time")
                elif mat["Data"][k].ndim == 2:
                    if "AHRSRotationMatrix" in k:
                        coords = {"dimRM": np.arange(9)}
                        dsi["AHRSRotationMatrix"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimRM"]
                        )
                    if "Magnetometer" in k:
                        coords = {"dimM": np.arange(3)}
                        dsi["Magnetometer"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimM"]
                        )
                    if "Accelerometer" in k:
                        coords = {"dimA": np.arange(3)}
                        dsi["Accelerometer"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimA"]
                        )
                    # only checks to see if cells match on first sample
                    elif dsi.attrs["data_type"] == "IBurst":
                        if mat["Data"][k].shape[1] == mat["Data"]["IBurst_NCells"][0]:
                            dsi[k.split("_")[1]] = xr.DataArray(
                                mat["Data"][k], dims=["time", "bindist"]
                            )
                    elif dsi.attrs["data_type"] == "IBurstHR":
                        if mat["Data"][k].shape[1] == mat["Data"]["IBurstHR_NCells"][0]:
                            dsi[k.split("_")[1]] = xr.DataArray(
                                mat["Data"][k], dims=["time", "bindist"]
                            )
                else:
                    print("still need to process", k, mat["Data"][k].shape)

        elif re.match("^Burst_", k) or re.match("^BurstHR_", k):
            if "_Time" not in k:
                if mat["Data"][k].ndim == 1:
                    dsb[k.split("_")[1]] = xr.DataArray(mat["Data"][k], dims="time")
                elif mat["Data"][k].ndim == 2:
                    if "AHRSRotationMatrix" in k:
                        coords = {"dimRM": np.arange(9)}
                        dsb["AHRSRotationMatrix"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimRM"]
                        )
                    if "Magnetometer" in k:
                        coords = {"dimM": np.arange(3)}
                        dsb["Magnetometer"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimM"]
                        )
                    if "Accelerometer" in k:
                        coords = {"dimA": np.arange(3)}
                        dsb["Accelerometer"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimA"]
                        )
                    # only checks to see if cells match on first sample
                    elif dsb.attrs["data_type"] == "Burst":
                        if mat["Data"][k].shape[1] == mat["Data"]["Burst_NCells"][0]:
                            dsb[k.split("_")[1]] = xr.DataArray(
                                mat["Data"][k], dims=["time", "bindist"]
                            )
                    elif dsb.attrs["data_type"] == "BurstHR":
                        if mat["Data"][k].shape[1] == mat["Data"]["BurstHR_NCells"][0]:
                            dsb[k.split("_")[1]] = xr.DataArray(
                                mat["Data"][k], dims=["time", "bindist"]
                            )
                else:
                    print("still need to process", k, mat["Data"][k].shape)

        elif re.match("^Echo1Bin1_", k):
            if "_Time" not in k:
                if mat["Data"][k].ndim == 1:
                    dse1[k.split("_")[2]] = xr.DataArray(mat["Data"][k], dims="time")
                elif mat["Data"][k].ndim == 2:
                    if "AHRSRotationMatrix" in k:
                        coords = {"dimRM": np.arange(9)}
                        dse1["AHRSRotationMatrix"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimRM"]
                        )
                    if "Magnetometer" in k:
                        coords = {"dimM": np.arange(3)}
                        dse1["Magnetometer"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimM"]
                        )
                    if "Accelerometer" in k:
                        coords = {"dimA": np.arange(3)}
                        dse1["Accelerometer"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimA"]
                        )
                    # only checks to see if cells match on first sample
                    elif (
                        mat["Data"][k].shape[1]
                        == mat["Data"][f"Echo1Bin1_{freq1}kHz_NCells"][0]
                    ):
                        dse1[k.split("_")[2]] = xr.DataArray(
                            mat["Data"][k], dims=["time", "bindist"]
                        )

                else:
                    print("still need to process", k, mat["Data"][k].shape)

        elif re.match("^Average_", k):
            if "_Time" not in k:
                if mat["Data"][k].ndim == 1:
                    dsa[k.split("_")[1]] = xr.DataArray(mat["Data"][k], dims="time")
                elif mat["Data"][k].ndim == 2:
                    if "AHRSRotationMatrix" in k:
                        coords = {"dimRM": np.arange(9)}
                        dsa["AHRSRotationMatrix"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimRM"]
                        )
                    if "Magnetometer" in k:
                        coords = {"dimM": np.arange(3)}
                        dsa["Magnetometer"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimM"]
                        )
                    if "Accelerometer" in k:
                        coords = {"dimA": np.arange(3)}
                        dsa["Accelerometer"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimA"]
                        )
                    # only checks to see if cells match on first sample
                    elif dsa.attrs["data_type"] == "Average":
                        if mat["Data"][k].shape[1] == mat["Data"]["Average_NCells"][0]:
                            dsa[k.split("_")[1]] = xr.DataArray(
                                mat["Data"][k], dims=["time", "bindist"]
                            )

                    else:
                        print("still need to process", k, mat["Data"][k].shape)

        elif re.match("^Alt_Average_", k):
            if "_Time" not in k:
                if mat["Data"][k].ndim == 1:
                    dsalt[k.split("_")[2]] = xr.DataArray(mat["Data"][k], dims="time")
                elif mat["Data"][k].ndim == 2:
                    if "AHRSRotationMatrix" in k:
                        coords = {"dimRM": np.arange(9)}
                        dsalt["AHRSRotationMatrix"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimRM"]
                        )
                    if "Magnetometer" in k:
                        coords = {"dimM": np.arange(3)}
                        dsalt["Magnetometer"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimM"]
                        )
                    if "Accelerometer" in k:
                        coords = {"dimA": np.arange(3)}
                        dsalt["Accelerometer"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimA"]
                        )
                    # only checks to see if cells match on first sample
                    elif dsalt.attrs["data_type"] == "Alt_Average":
                        if mat["Data"][k].shape[1] == mat["Data"]["Average_NCells"][0]:
                            dsalt[k.split("_")[2]] = xr.DataArray(
                                mat["Data"][k], dims=["time", "bindist"]
                            )

                    else:
                        print("still need to process", k, mat["Data"][k].shape)
        else:
            print("missing variable:", k)

    for ds in ds_dict:
        read_config_mat(mat, ds_dict[ds])

    for ds in ds_dict:
        add_descriptions(mat, ds_dict[ds])

    for ds in ds_dict:
        add_units(mat, ds_dict[ds])
        add_transmatrix(mat, ds_dict[ds])

    return ds_dict


def read_config_mat(mat, ds):
    for k in mat["Config"]:
        if re.search("_Beam2xyz$", k):
            ds.attrs[f"SIG{k}"] = str(mat["Config"][k])
        else:
            ds.attrs[f"SIG{k}"] = mat["Config"][k]


def add_descriptions(mat, ds):
    for k in mat["Descriptions"]:
        var = k.split("_")[1]
        if var in ds:
            if "long_name" not in ds[var].attrs:
                ds[var].attrs["long_name"] = mat["Descriptions"][k]


def add_units(mat, ds):
    for k in mat["Units"]:
        var = k.split("_")[1]
        if var in ds:
            if "units" not in ds[var].attrs:
                ds[var].attrs["units"] = mat["Units"][k]


def add_transmatrix(mat, ds):
    for k in mat["Config"]:
        if f"{ds.attrs['data_type']}_Beam2xyz" in k:
            var = k.split("_")[-1]
            ds[var] = xr.DataArray(mat["Config"][k])


def mat_to_cdf(metadata):
    """
    Load .mat files exported from Signature software and process to .cdf
    """

    tic = time.time()
    basefile = metadata["basefile"]

    if "prefix" in metadata:
        prefix = metadata["prefix"]
    else:
        prefix = ""
    basefile = prefix + basefile

    if "outdir" in metadata:
        outdir = metadata["outdir"]
    else:
        outdir = ""

    utils.check_valid_globalatts_metadata(metadata)
    aqdutils.check_valid_config_metadata(metadata)

    print("Loading from Matlab files; this may take a while for large datasets")

    def loadpar(f):
        dsd = load_mat_file(f)
        filstub = os.path.normpath(f).split("\\")[-1]
        bf = basefile.split("/")[-1]
        num = filstub.split(f"{bf}_")[-1].split(".mat")[0]
        for dsn in dsd:
            ds = dsd[dsn]
            ds = utils.write_metadata(ds, metadata)
            ds = utils.ensure_cf(ds)
            cdf_filename = (
                prefix
                + outdir
                + ds.attrs["filename"]
                + "-"
                + ds.attrs["data_type"]
                + "-"
                + num
                + "-raw.cdf"
            )
            ds.to_netcdf(cdf_filename)
            print(f"Finished writing data to {cdf_filename}")

    matfiles = glob.glob(f"{basefile}_*.mat")
    # print(matfiles)
    if len(matfiles) > 1:
        Parallel(n_jobs=-1, verbose=10)(delayed(loadpar)(f) for f in matfiles)
    else:
        loadpar(matfiles[0])

    dsd = load_mat_file(matfiles[0])  # get minimal dsd file for below
    cdffiles = glob.glob(f"{prefix}{outdir}*-raw.cdf")
    ds = xr.open_dataset(cdffiles[0])  # get minimal ds for below

    # read in Burst -raw.cdf and make one combined per data_type
    if "chunks" in ds.attrs:
        chunksizes = dict(zip(ds.attrs["chunks"][::2], ds.attrs["chunks"][1::2]))
        for key in chunksizes:
            if isinstance(chunksizes[key], str):
                chunksizes[key] = int(chunksizes[key])
        print(f"Using user specified chunksizes = {chunksizes}")
    else:
        chunksizes = {"time": 200000, "bindist": 48}
        print(f"Using default chunksizes = {chunksizes}")

    for k in dsd:
        # dsb = dsd["dsb"]
        fin = outdir + f"*-{dsd[k].attrs['data_type']}-*.cdf"
        print(k)
        try:
            ds = xr.open_mfdataset(fin, parallel=True, chunks=chunksizes)
            ds = aqdutils.check_attrs(ds, inst_type="SIG")
            if "Beam2xyz" in ds:
                if "time" in ds["Beam2xyz"].dims:
                    ds["Beam2xyz"] = ds["Beam2xyz"].isel(time=0, drop=True)
            # write out all into single -raw.cdf files per data_type
            if (
                ds.attrs["data_type"].lower() == "burst"
                or ds.attrs["data_type"].lower() == "bursthr"
            ):
                ftype = "burst"
            elif (
                ds.attrs["data_type"].lower() == "iburst"
                or ds.attrs["data_type"].lower() == "ibursthr"
            ):
                ftype = "iburst"
            elif ds.attrs["data_type"].lower() == "echosounder":
                ftype = "echo1"
            elif ds.attrs["data_type"].lower() == "burstrawaltimeter":
                ftype = "burstrawalt"
            elif ds.attrs["data_type"].lower() == "average":
                ftype = "avgd"
            elif ds.attrs["data_type"].lower() == "alt_average":
                ftype = "altavgd"

            cdf_filename = prefix + ds.attrs["filename"] + f"_{ftype}-raw.cdf"
            print(f"writing {ftype} to netcdf")
            delayed_obj = ds.to_netcdf(cdf_filename, compute=False)
            with ProgressBar():
                delayed_obj.compute()

            print(f"Finished writing data to {cdf_filename}")

        except ValueError as ve:
            print(
                f"Failed to read or write netCDF file(s) for data_type {dsd[k].attrs['data_type']} to make composite -raw.cdf file, due to Value Error: '{ve}'"
            )

    toc = time.time()
    etime = round(toc - tic, 0)
    print(f"elapsed time = {etime}")
