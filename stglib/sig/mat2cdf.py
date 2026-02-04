import datetime as dt
import glob
import re
import time

import numpy as np
import pandas as pd
import xarray as xr
from dask.diagnostics import ProgressBar
from tqdm import tqdm

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

    if "Alt_Plan_BurstEnabled" in mat["Config"]:
        if mat["Config"]["Alt_Plan_BurstEnabled"] == "True":
            bindistAltBurst = (
                mat["Config"]["Alt_Burst_BlankingDistance"]
                + mat["Config"]["Alt_Burst_CellSize"] / 2
                + mat["Config"]["Alt_Burst_CellSize"]
                * np.arange(mat["Config"]["Alt_Burst_NCells"])
            )

        if (
            mat["Config"]["Alt_Burst_HighResolution"] == "True"
            or mat["Config"]["Alt_Burst_HighResolution5"] == "True"
        ):
            bindistAltBurstHR = (
                mat["Config"]["Alt_BurstHR_BlankingDistance"]
                + mat["Config"]["Alt_BurstHR_CellSize"] / 2
                + mat["Config"]["Alt_BurstHR_CellSize"]
                * np.arange(mat["Config"]["Alt_BurstHR_NCells"])
            )

        if mat["Config"]["Alt_Burst_EchoSounder"] == "True":
            bindistAltECHO = (
                mat["Config"]["Alt_EchoSounder_BlankingDistance"]
                + mat["Config"]["Alt_EchoSounder_CellSize"] / 2
                + mat["Config"]["Alt_EchoSounder_CellSize"]
                * np.arange(mat["Config"]["Alt_EchoSounder_NCells"])
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

    # Alt Average
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

    # Alt Burst
    if "Alt_Plan_BurstEnabled" in mat["Config"]:
        if (
            mat["Config"]["Alt_Plan_BurstEnabled"] == "True"
            and mat["Config"]["Alt_Burst_RawAltimeter"] == 1
            and mat["Config"]["Alt_Burst_Altimeter"] == "True"
        ):
            # Alt BurstRawAltimeter
            if "Alt_BurstRawAltimeter_Time" in mat["Data"]:
                dsaltbra = xr.Dataset()
                dsaltbra["time"] = xr.DataArray(
                    [
                        matlab2datetime(x)
                        for x in mat["Data"]["Alt_BurstRawAltimeter_Time"]
                    ],
                    dims="time",
                )
                dsaltbra["time"] = pd.DatetimeIndex(dsaltbra["time"])
                dsaltbra["time"] = pd.DatetimeIndex(dsaltbra["time"])
                dsaltbra.attrs["data_type"] = "Alt_BurstRawAltimeter"
                ds_dict["dsaltbra"] = dsaltbra

        if (
            mat["Config"]["Alt_Plan_BurstEnabled"] == "True"
            and mat["Config"]["Alt_Burst_NBeams"] == 5
        ):
            # Alt IBurst
            if (
                mat["Config"]["Alt_Burst_HighResolution5"] == "True"
                and "Alt_IBurstHR_Time" in mat["Data"]
            ):
                dsaltbi = xr.Dataset()
                dsaltbi["time"] = pd.DatetimeIndex(
                    xr.DataArray(
                        [matlab2datetime(x) for x in mat["Data"]["Alt_IBurstHR_Time"]],
                        dims="time",
                    )
                )
                dsaltbi["time"] = pd.DatetimeIndex(dsaltbi["time"])
                dsaltbi["bindist"] = xr.DataArray(bindistAltBurstHR, dims="bindist")
                dsaltbi.attrs["data_type"] = "Alt_IBurstHR"
                ds_dict["dsaltbi"] = dsaltbi

            elif "Alt_IBurst_Time" in mat["Data"]:
                dsaltbi = xr.Dataset()
                dsaltbi["time"] = pd.DatetimeIndex(
                    xr.DataArray(
                        [matlab2datetime(x) for x in mat["Data"]["Alt_IBurst_Time"]],
                        dims="time",
                    )
                )
                dsaltbi["time"] = pd.DatetimeIndex(dsaltbi["time"])
                dsaltbi["bindist"] = xr.DataArray(bindistAltBurst, dims="bindist")
                dsaltbi.attrs["data_type"] = "Alt_IBurst"
                ds_dict["dsaltbi"] = dsaltbi

        if mat["Config"]["Alt_Plan_BurstEnabled"] == "True":

            # Alt Burst
            if (
                mat["Config"]["Alt_Burst_HighResolution"] == "True"
                and "Alt_BurstHR_Time" in mat["Data"]
            ):
                dsaltb = xr.Dataset()
                dsaltb["time"] = pd.DatetimeIndex(
                    xr.DataArray(
                        [matlab2datetime(x) for x in mat["Data"]["Alt_BurstHR_Time"]],
                        dims="time",
                    )
                )
                dsaltb["time"] = pd.DatetimeIndex(dsaltb["time"])
                dsaltb["bindist"] = xr.DataArray(bindistAltBurstHR, dims="bindist")
                dsaltb.attrs["data_type"] = "Alt_BurstHR"
                ds_dict["dsaltb"] = dsaltb

            elif "Alt_Burst_Time" in mat["Data"]:
                dsaltb = xr.Dataset()
                dsaltb["time"] = pd.DatetimeIndex(
                    xr.DataArray(
                        [matlab2datetime(x) for x in mat["Data"]["Alt_Burst_Time"]],
                        dims="time",
                    )
                )
                dsaltb["time"] = pd.DatetimeIndex(dsaltb["time"])
                dsaltb["bindist"] = xr.DataArray(bindistAltBurst, dims="bindist")
                dsaltb.attrs["data_type"] = "Alt_Burst"
                ds_dict["dsaltb"] = dsaltb

        if (
            mat["Config"]["Alt_Plan_BurstEnabled"] == "True"
            and mat["Config"]["Alt_Burst_EchoSounder"] == "True"
        ):
            # echo1 data - only handling echo1 data to start
            if "Alt_EchoSounder_Frequency1" in mat["Config"]:
                freq1 = mat["Config"]["Alt_EchoSounder_Frequency1"]
                if f"Alt_Echo1Bin1_{freq1}kHz_Time" in mat["Data"]:
                    dsalte1 = xr.Dataset()
                    dsalte1["time"] = pd.DatetimeIndex(
                        xr.DataArray(
                            [
                                matlab2datetime(x)
                                for x in mat["Data"][f"Alt_Echo1Bin1_{freq1}kHz_Time"]
                            ],
                            dims="time",
                        )
                    )
                    dsalte1["time"] = pd.DatetimeIndex(dsalte1["time"])
                    dsalte1["bindist"] = xr.DataArray(bindistAltECHO, dims="bindist")
                    dsalte1.attrs["data_type"] = "Alt_EchoSounder"
                    ds_dict["dsalte1"] = dsalte1

    # add data to data sets
    for k in mat["Data"]:
        if re.match("^BurstRawAltimeter_", k):
            if "_Time" not in k:
                if mat["Data"][k].ndim == 1:
                    dsbra[k.split("_")[1]] = xr.DataArray(mat["Data"][k], dims="time")
                elif mat["Data"][k].ndim == 2:
                    # print("still need to process", k, mat["Data"][k].shape)
                    pass
        elif re.match("^IBurst_", k) or re.match("^IBurstHR_", k):
            if "_Time" not in k:
                if mat["Data"][k].ndim == 1:
                    dsi[k.split("_")[1]] = xr.DataArray(mat["Data"][k], dims="time")
                elif mat["Data"][k].ndim == 2:
                    if "AHRSRotationMatrix" in k:
                        dsi["AHRSRotationMatrix"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimRM"]
                        )
                    if "Magnetometer" in k:
                        dsi["Magnetometer"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimM"]
                        )
                    if "Accelerometer" in k:
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
                        dsb["AHRSRotationMatrix"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimRM"]
                        )
                    if "Magnetometer" in k:
                        dsb["Magnetometer"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimM"]
                        )
                    if "Accelerometer" in k:
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
                        dse1["AHRSRotationMatrix"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimRM"]
                        )
                    if "Magnetometer" in k:
                        dse1["Magnetometer"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimM"]
                        )
                    if "Accelerometer" in k:
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
                        dsa["AHRSRotationMatrix"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimRM"]
                        )
                    if "Magnetometer" in k:
                        dsa["Magnetometer"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimM"]
                        )
                    if "Accelerometer" in k:
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
                        dsalt["AHRSRotationMatrix"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimRM"]
                        )
                    if "Magnetometer" in k:
                        dsalt["Magnetometer"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimM"]
                        )
                    if "Accelerometer" in k:
                        dsalt["Accelerometer"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimA"]
                        )
                    # only checks to see if cells match on first sample
                    elif dsalt.attrs["data_type"] == "Alt_Average":
                        if (
                            mat["Data"][k].shape[1]
                            == mat["Data"]["Alt_Average_NCells"][0]
                        ):
                            dsalt[k.split("_")[2]] = xr.DataArray(
                                mat["Data"][k], dims=["time", "bindist"]
                            )

                    else:
                        print("still need to process", k, mat["Data"][k].shape)

        elif re.match("^Alt_BurstRawAltimeter_", k):
            if "_Time" not in k:
                if mat["Data"][k].ndim == 1:
                    dsaltbra[k.split("_")[2]] = xr.DataArray(
                        mat["Data"][k], dims="time"
                    )
                elif mat["Data"][k].ndim == 2:
                    if "AHRSRotationMatrix" in k:
                        dsaltbra["AHRSRotationMatrix"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimRM"]
                        )
                    if "Magnetometer" in k:
                        dsaltbra["Magnetometer"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimM"]
                        )
                    if "Accelerometer" in k:
                        dsaltbra["Accelerometer"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimA"]
                        )
                    # print("still need to process", k, mat["Data"][k].shape)
                    # pass
        elif re.match("^Alt_IBurst_", k) or re.match("^Alt_IBurstHR_", k):
            if "_Time" not in k:
                if mat["Data"][k].ndim == 1:
                    dsaltbi[k.split("_")[2]] = xr.DataArray(mat["Data"][k], dims="time")
                elif mat["Data"][k].ndim == 2:
                    if "AHRSRotationMatrix" in k:
                        dsaltbi["AHRSRotationMatrix"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimRM"]
                        )
                    if "Magnetometer" in k:
                        dsaltbi["Magnetometer"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimM"]
                        )
                    if "Accelerometer" in k:
                        dsaltbi["Accelerometer"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimA"]
                        )
                    # only checks to see if cells match on first sample
                    elif dsaltbi.attrs["data_type"] == "Alt_IBurst":
                        if (
                            mat["Data"][k].shape[1]
                            == mat["Data"]["Alt_IBurst_NCells"][0]
                        ):
                            dsaltbi[k.split("_")[2]] = xr.DataArray(
                                mat["Data"][k], dims=["time", "bindist"]
                            )
                    elif dsaltbi.attrs["data_type"] == "Alt_IBurstHR":
                        if (
                            mat["Data"][k].shape[1]
                            == mat["Data"]["Alt_IBurstHR_NCells"][0]
                        ):
                            dsaltbi[k.split("_")[2]] = xr.DataArray(
                                mat["Data"][k], dims=["time", "bindist"]
                            )
                else:
                    print("still need to process", k, mat["Data"][k].shape)

        elif re.match("^Alt_Burst_", k) or re.match("^Alt_BurstHR_", k):
            if "_Time" not in k:
                if mat["Data"][k].ndim == 1:
                    dsaltb[k.split("_")[2]] = xr.DataArray(mat["Data"][k], dims="time")
                elif mat["Data"][k].ndim == 2:
                    if "AHRSRotationMatrix" in k:
                        dsaltb["AHRSRotationMatrix"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimRM"]
                        )
                    if "Magnetometer" in k:
                        dsaltb["Magnetometer"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimM"]
                        )
                    if "Accelerometer" in k:
                        dsaltb["Accelerometer"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimA"]
                        )
                    # only checks to see if cells match on first sample
                    elif dsaltb.attrs["data_type"] == "Alt_Burst":
                        if (
                            mat["Data"][k].shape[1]
                            == mat["Data"]["Alt_Burst_NCells"][0]
                        ):
                            dsaltb[k.split("_")[2]] = xr.DataArray(
                                mat["Data"][k], dims=["time", "bindist"]
                            )
                    elif dsaltb.attrs["data_type"] == "Alt_BurstHR":
                        if (
                            mat["Data"][k].shape[1]
                            == mat["Data"]["Alt_BurstHR_NCells"][0]
                        ):
                            dsaltb[k.split("_")[2]] = xr.DataArray(
                                mat["Data"][k], dims=["time", "bindist"]
                            )
                else:
                    print("still need to process", k, mat["Data"][k].shape)

        elif re.match("^Alt_Echo1Bin1_", k):
            if "_Time" not in k:
                if mat["Data"][k].ndim == 1:
                    dsalte1[k.split("_")[2]] = xr.DataArray(mat["Data"][k], dims="time")
                elif mat["Data"][k].ndim == 2:
                    if "AHRSRotationMatrix" in k:
                        dsalte1["AHRSRotationMatrix"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimRM"]
                        )
                    if "Magnetometer" in k:
                        dsalte1["Magnetometer"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimM"]
                        )
                    if "Accelerometer" in k:
                        dsalte1["Accelerometer"] = xr.DataArray(
                            mat["Data"][k], dims=["time", "dimA"]
                        )
                    # only checks to see if cells match on first sample
                    elif (
                        mat["Data"][k].shape[1]
                        == mat["Data"][f"Alt_Echo1Bin1_{freq1}kHz_NCells"][0]
                    ):
                        dsalte1[k.split("_")[2]] = xr.DataArray(
                            mat["Data"][k], dims=["time", "bindist"]
                        )

                else:
                    print("still need to process", k, mat["Data"][k].shape)

        else:
            print("missing variable:", k)

    for ds in ds_dict:
        # if burst plan data find sample_mode ['BURST' or 'CONTINUOUS']
        if ds == "dsi" or ds == "dsb" or ds == "dse1":
            if (
                mat["Config"]["Plan_BurstInterval"]
                * mat["Config"]["Burst_SamplingRate"]
                == mat["Config"]["Burst_NSample"]
            ):
                ds_dict[ds].attrs["sample_mode"] = "CONTINUOUS"
            else:
                ds_dict[ds].attrs["sample_mode"] = "BURST"

        elif ds == "dsa" or ds == "dsalt":
            ds_dict[ds].attrs["sample_mode"] = "AVERAGE"

        elif ds == "dsaltbi" or ds == "dsaltb" or ds == "dsalte1":
            if (
                mat["Config"]["Alt_Plan_BurstInterval"]
                * mat["Config"]["Alt_Burst_SamplingRate"]
                == mat["Config"]["Alt_Burst_NSample"]
            ):
                ds_dict[ds].attrs["sample_mode"] = "CONTINUOUS"
            else:
                ds_dict[ds].attrs["sample_mode"] = "BURST"

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

    utils.check_valid_globalatts_metadata(metadata)
    aqdutils.check_valid_config_metadata(metadata)

    print("Loading from Matlab files; this may take a while for large datasets")

    def make_sig_dsd(basefile):
        matfiles = glob.glob(f"{basefile}_*.mat")
        dsd = load_mat_file(matfiles[0])  # get minimal dsd file for below
        ds_dict = {}
        for k in dsd:
            ds_dict[k] = []

        for f in tqdm(matfiles):
            dsd = load_mat_file(f)

            for dsn in dsd:
                ds_dict[dsn].append(dsd[dsn].chunk({"time": 300000}))

        return ds_dict

    dsd = make_sig_dsd(basefile)

    for k in dsd:
        # dsb = dsd["dsb"]
        ds = xr.concat(dsd[k], dim="time").chunk({"time": 300000})
        ds = utils.write_metadata(ds, metadata)
        ds = utils.ensure_cf(ds)

        ds = aqdutils.check_attrs(ds, inst_type="SIG")

        print(
            f"sorting {k} dataset by time, this can take some time to complete for very large datasets"
        )
        ds = ds.sortby("time")
        print(f"sorting by time completed for {k} dataset")

        if "Beam2xyz" in ds:
            if "time" in ds["Beam2xyz"].dims:
                ds["Beam2xyz"] = ds["Beam2xyz"].isel(time=0, drop=True)
        # write out all into single -raw.cdf files per data_type
        if ds.attrs["data_type"] in ["Burst", "BurstHR"]:
            ftype = "burst"
        elif ds.attrs["data_type"] in ["IBurst", "IBurstHR"]:
            ftype = "iburst"
        elif ds.attrs["data_type"] == "EchoSounder":
            ftype = "echo1"
        elif ds.attrs["data_type"] == "BurstRawAltimeter":
            ftype = "burstrawalt"
        elif ds.attrs["data_type"] == "Average":
            ftype = "avgd"
        elif ds.attrs["data_type"] == "Alt_Average":
            ftype = "altavgd"
        elif ds.attrs["data_type"] in ["Alt_Burst", "Alt_BurstHR"]:
            ftype = "altburst"
        elif ds.attrs["data_type"] in ["Alt_IBurst", "Alt_IBurstHR"]:
            ftype = "altiburst"
        elif ds.attrs["data_type"] == "Alt_EchoSounder":
            ftype = "altecho1"

        elif ds.attrs["data_type"] == "Alt_BurstRawAltimeter":
            ftype = "altburstrawalt"

        cdf_filename = prefix + ds.attrs["filename"] + f"_{ftype}-raw.cdf"
        print(f"writing {ftype} to netcdf")
        delayed_obj = ds.to_netcdf(cdf_filename, compute=False)
        with ProgressBar():
            delayed_obj.compute()

        print(f"Finished writing data to {cdf_filename}")

    toc = time.time()
    etime = round(toc - tic, 0)
    print(f"elapsed time = {etime}")
