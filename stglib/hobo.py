from __future__ import division, print_function
import pandas as pd
import xarray as xr


def read_hobo(filnam, skiprows=1, skipfooter=0):
    hobo =  pd.read_csv('/Volumes/Backstaff/field/1104/hobo3353atmos/3353-1104-atmos.csv',
                      usecols=[0, 1, 2, 3],
                      names=['#','datetime','abspres_kPa','temp_C'],
                      engine='python',
                      skiprows=skiprows,
                      skipfooter=skipfooter)
    hobo['time'] = pd.to_datetime(hobo['datetime'])
    hobo['abspres_dbar'] = hobo['abspres_kPa']/10
    hobo.set_index('time', inplace=True)

    return xr.Dataset(hobo)
