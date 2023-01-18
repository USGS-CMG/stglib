import sys

import numpy as np
import pandas as pd
import scipy.io
import xarray as xr

from stglib.core.utils import loadmat

sys.path.append("/Users/dnowacki/Documents/python")
import matlabtools

# %%
smat = scipy.io.loadmat(
    "/Users/dnowacki/Downloads/sig1k/S100151A012_CSF20CHT_1.mat", simplify_cells=True
)
mat = loadmat("/Users/dnowacki/Downloads/sig1k/S100151A012_CSF20CHT_1.mat")
mat["__version__"]
# %%
smat["Config"].keys()
mat["Config"].keys()
np.array_equal(smat["Config"], mat["Config"])
np.all(smat["Config"] == mat["Config"])
np.all(mat["Config"] == smat["Config"])
# %%

# %%

bindist = (
    mat["Config"]["Burst_BlankingDistance"]
    + mat["Config"]["Burst_CellSize"] / 2
    + mat["Config"]["Burst_CellSize"] * np.arange(mat["Config"]["Burst_NCells"])
)

for k in mat["Config"]:
    ds.attrs[k] = mat["Config"][k]
# %%
dsi = xr.Dataset()
dsi["time"] = pd.DatetimeIndex(
    xr.DataArray(
        [matlabtools.matlab2datetime(x) for x in mat["Data"]["IBurst_Time"]],
        dims="time",
    )
)
dsi["time"] = pd.DatetimeIndex(dsi["time"])
ds = xr.Dataset()
ds["time"] = xr.DataArray(
    [matlabtools.matlab2datetime(x) for x in mat["Data"]["BurstRawAltimeter_Time"]],
    dims="time",
)
ds["time"] = pd.DatetimeIndex(ds["time"])
ds["time"] = pd.DatetimeIndex(ds["time"])
for k in mat["Data"]:
    # if 'BurstRawAltimeter' in k:
    #     if '_Time' not in k:
    #         if mat['Data'][k].ndim == 1:
    #             ds[k.split('_')[1]] = xr.DataArray(mat['Data'][k], dims='time')
    #         elif mat['Data'][k].ndim == 2:
    #             print(k)
    #             print(mat['Data'][k].shape)
    if "IBurst" in k:
        if "_Time" not in k:
            if mat["Data"][k].ndim == 1:
                dsi[k.split("_")[1]] = xr.DataArray(mat["Data"][k], dims="time")
            elif mat["Data"][k].ndim == 2:
                # only checks to see if cells match on first sample
                if mat["Data"][k].shape[1] == mat["Data"]["IBurst_NCells"][0]:
                    dsi[k.split("_")[1]] = xr.DataArray(
                        mat["Data"][k], dims=["time", "cells"]
                    )
                    print(k, mat["Data"][k].shape)
# %%
dsi.time.diff(dim="time") / pd.Timedelta("1s")
dsi["VelBeam5"].plot()
# %%
print(dsi.NBeams)
for k in dsi:
    print(k)
# %%
mat["Data"]
mat["Units"]
# %%
ds.attrs
# %%
for k in mat.keys():
    print("***", k)
    if type(mat[k]) == dict:
        print(mat[k].keys())
    else:
        print(mat[k])
