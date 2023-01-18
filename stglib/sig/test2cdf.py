import sys

import numpy as np

sys.path.append("/Users/dnowacki/Documents/python")
import matlabtools
import xarray as xr

# %%
fil = "/Volumes/Backstaff/scratch/CSF20CHT01sig/rawdata/S100151A012_CSF20CHT/S100151A012_CSF20CHT_1.mat"
# %%


# %%
mat = matlabtools.loadmat(fil)
# %%
mat["Descriptions"]
(
    mat["Config"]["Burst_BlankingDistance"]
    + mat["Config"]["Burst_CellSize"] / 2
    + mat["Config"]["Burst_CellSize"] * np.arange(mat["Config"]["Burst_NCells"])
)

mat["Config"]
