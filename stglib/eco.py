from __future__ import division, print_function
import pandas as pd
import xarray as xr


def read_par(filnam, spb=False, skiprows=None, skipfooter=0):
    """Read data from a WET Labs PAR csv file into an xarray
    Dataset.

    Parameters
    ----------
    filnam : string
        The filename
    spb: bool, optional
        Samples per burst if using burst sampling
    skiprows : int, optional
        How many header rows to skip. Default None
    skipfooter : into, optional
        How many footer rows to skip. Default 0

    Returns
    -------
    xarray.Dataset
        An xarray Dataset of the PAR data
    """

    names = ['date', 'time', 'counts']

    par = read_eco_csv(filnam, names, skiprows=skiprows, skipfooter=skipfooter)

    return eco_pd_to_xr(par, spb=spb)


def read_ntu(filnam, spb=False, skiprows=None, skipfooter=0):
    """Read data from a WET Labs NTU csv file into an xarray
    Dataset.

    Parameters
    ----------
    filnam : string
        The filename
    spb: bool, optional
        Samples per burst if using burst sampling
    skiprows : int, optional
        How many header rows to skip. Default None
    skipfooter : into, optional
        How many footer rows to skip. Default 0

    Returns
    -------
    xarray.Dataset
        An xarray Dataset of the PAR data
    """

    names = ['date', 'time', 'a', 'counts', 'b']

    ntu = read_eco_csv(filnam, names, skiprows=skiprows, skipfooter=skipfooter)

    return eco_pd_to_xr(ntu, spb=spb)


def read_eco_csv(filnam, names, skiprows=None, skipfooter=0):

    return pd.read_csv(filnam,
                      sep='\t',
                      names=names,
                      parse_dates=[['date', 'time']],
                      infer_datetime_format=True,
                      engine='python',
                      skiprows=skiprows,
                      skipfooter=skipfooter)


def eco_pd_to_xr(df, spb=False):

    if spb:
        times = df['date_time'].values.reshape((-1, spb))[:, int(spb/2)] # get middle time
        counts = df['counts'].values.reshape((-1, spb))
        sample = range(spb)

        ds = xr.Dataset({'time': ('time', times),
                         'counts': (['time', 'sample'], counts),
                         'sample': ('sample', sample)})
    else:
        times = df['date_time']
        counts = df['counts']

        ds = xr.Dataset({'time': ('time', times),
                         'counts': ('time', counts)})

    return ds
