from __future__ import division, print_function
import pandas as pd
import xarray as xr


def read_hobo(filnam, skiprows=1, skipfooter=0):
    """Read data from an Onset HOBO pressure sensor .csv file into an xarray
    Dataset.

    Parameters
    ----------
    filnam : string
        The filename
    skiprows : int, optional
        How many header rows to skip. Default 1
    skipfooter : int, optional
        How many footer rows to skip. Default 0

    Returns
    -------
    xarray.Dataset
        An xarray Dataset of the HOBO data
    """
    hobo =  pd.read_csv(filnam,
                      usecols=[0, 1, 2, 3],
                      names=['#','datetime','abspres_kPa','temp_C'],
                      engine='python',
                      skiprows=skiprows,
                      skipfooter=skipfooter)
    hobo['time'] = pd.to_datetime(hobo['datetime'])
    hobo['abspres_dbar'] = hobo['abspres_kPa']/10
    hobo.set_index('time', inplace=True)

    return xr.Dataset(hobo)
