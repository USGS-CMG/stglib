import dolfyn as dlfn
from dolfyn.adv import api
import pandas as pd

# takes almost 6 minutes to read 400MB .VEC file
fildir = "/Volumes/Backstaff/scratch/WTW21MES/WTW21MES01/"
# %time ds = dlfn.read(f'{fildir}WTWMES02.VEC')
ds = dlfn.read(f"{fildir}WTWMES02.VEC")
# %%
# ds_orig = dlfn.read('/Users/dnowacki/projects/dolfyn/dolfyn/example_data/burst_mode01.VEC')
# ds = dlfn.read('/Users/dnowacki/projects/dolfyn/dolfyn/example_data/burst_mode01.VEC')
# dlfn.rotate2(ds, 'earth', inplace=True)

# %%

# First set the magnetic declination
dlfn.set_declination(ds, declin=10, inplace=True)  # declination points 10 degrees East

# Rotate that data from the instrument to earth frame (ENU):
dlfn.rotate2(ds, "earth", inplace=True)
# %%
# ds['time'] = ds.time - (ds.time[0] - ds.time[0].dt.round('1s')) # 1 msec offset
dlfn.save(ds, "/Volumes/Backstaff/scratch/WTW21MES/WTW21MES01/dolfyn_raw.nc")
ds.to_netcdf("/Volumes/Backstaff/scratch/WTW21MES/WTW21MES01/dolfyn_raw2.nc")
# %%

# %%

binner = api.ADVBinner(n_bin=ds.attrs["DutyCycle_NBurst"], fs=ds.fs)
ds_binned = binner(ds, freq_units="Hz")
# %%
for var in ["tke_vec", "stress", "psd"]:
    if var in ds_binned:
        ds_binned = ds_binned.drop_vars(var)

for k in ds_binned:
    if "dir" in ds_binned[k].coords:
        for d in ds_binned.dir:
            ds_binned[f"{k}_{d.values}"] = ds_binned[k].sel(dir=d)
            ds_binned[f"{k}_{d.values}"].attrs = ds_binned[k].sel(dir=d).attrs
        ds_binned = ds_binned.drop(k)

usedcoords = []
for k in ds_binned.data_vars:
    for c in ds_binned[k].coords:
        usedcoords.append(c)

usedcoords = set(usedcoords)

for k in ds_binned.coords:
    if k not in usedcoords:
        ds_binned = ds_binned.drop(k)

standard_names = {
    "c_sound": "speed_of_sound_in_sea_water",
    "heading": "platform_orientation",
    "pitch": "platform_pitch",
    "roll": "platform_roll",
    "temp": "sea_water_temperature",
    "vel_E": "eastward_sea_water_velocity",
    "vel_N": "northward_sea_water_velocity",
    "vel_U": "upward_sea_water_velocity",
    "pressure": "sea_water_pressure",
}

long_names = {
    "U_std": "Velocity standard deviation",
    "batt": "Battery voltage",
    "error": "",
    "status": "",
    "amp_E": "Beam 1 amplitude",
    "amp_N": "Beam 2 amplitude",
    "amp_U": "Beam 3 amplitude",
    "corr_E": "Beam 1 correlation",
    "corr_N": "Beam 2 correlation",
    "corr_U": "Beam 3 correlation",
    "orientation_down": "",
    "orientmat": "",
}

for k in standard_names:
    if standard_names[k] != "":
        ds_binned[k].attrs["standard_name"] = standard_names[k]

for k in long_names:
    if long_names[k] != "":
        ds_binned[k].attrs["long_name"] = long_names[k]

for var in ["heading", "pitch", "roll"]:
    ds_binned[var].attrs["units"] = "degree"
ds_binned["temp"].attrs["units"] = "degree_C"
import numpy as np

np.unique(ds.time.dt.round("1ms").diff(dim="time"))
np.unique(ds.time.diff(dim="time"))
ds_binned.time[:5]
ds_binned.time[-5:]
# %%
dlfn.save(ds_binned, f"{fildir}dolfyn_binned.nc")
ds_binned.to_netcdf(f"{fildir}dolfyn_binned.cdf")
