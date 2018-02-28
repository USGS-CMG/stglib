from __future__ import division, print_function
import pandas as pd
import xarray as xr

def read_exo(filnam, skiprows=25):
    """Read data from a YSI EXO multiparameter sonde .xlsx file into an xarray
    Dataset.

    Parameters
    ----------
    filnam : string
        The filename
    skiprows : int, optional
        How many header rows to skip. Default 25

    Returns
    -------
    xarray.Dataset
        An xarray Dataset of the EXO data
    """

    exo = pd.read_excel(filnam,
                        skiprows=skiprows,
                        infer_datetime_format=True,
                        parse_dates=[['Date (MM/DD/YYYY)', 'Time (HH:MM:SS)']])
    exo.rename(columns={'Date (MM/DD/YYYY)_Time (HH:MM:SS)': 'time'}, inplace=True)
    exo.set_index('time', inplace=True)

    return xr.Dataset(exo)
