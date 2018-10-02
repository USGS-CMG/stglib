import warnings
import pandas as pd
import xarray as xr
import numpy as np
import scipy.signal
from .core import utils
from .aqd import qaqc

def read_argonaut(filbase):
    return pd.read_csv(filnam,
                       skiprows=skiprows,
                       infer_datetime_format=True,
                       parse_dates=['Date and Time'],
                       encoding=encoding)

def read_dat_raw(filnam):
    df = pd.read_csv(filnam,
                     delim_whitespace=True,
                     parse_dates={'time': ['Year',
                                           'Month',
                                           'Day',
                                           'Hour',
                                           'Minute',
                                           'Second']},
                     date_parser=qaqc.date_parser)
    df.set_index('time', inplace=True)
    return df

def read_dat(filnam):
    return read_dat_raw(filnam)

def read_vel(filnam):
    df = pd.read_csv(filnam,
                     delim_whitespace=True,
                     header=[0,1],
                     parse_dates=[[1,2,3,4,5,6]],
                     date_parser=qaqc.date_parser)
    df.rename(columns='_'.join, inplace=True)
    df.rename(columns={"(_'_Y_'_,_ _'_(_)_'_)___(_'_M_'_,_ _'_(_)_'_)___(_'_D_'_,_ _'_(_)_'_)___(_'_H_'_,_ _'_(_)_'_)___(_'_M_'_,_ _'_(_)_'_)_._1___(_'_S_'_,_ _'_(_)_'_)":
            'time'},
            index=str,
            inplace=True)
    # df.columns = df.columns.str.replace('(cm/s)', '')
    # df.columns = df.columns.str.replace('(deg)', '')
    # df.columns = df.columns.str.replace('()', '')
    df.columns = df.columns.str.replace(r"\(.*\)","")
    df.set_index('time',inplace=True)
    ds = xr.Dataset(df)
    
    return df

def read_ctl(filnam):
    with open( filnam) as f:
        row = ''
        while 'Argonaut ASCII data file Long format is as follows' not in row:
            row = f.readline().rstrip()

        while 'Flow data file format is as follows' not in row:
            row = f.readline().rstrip()

        f.close()

# parse_dates={'datetime': [2, 0, 1, 3, 4, 5]},
# date_parser=qaqc.date_parser,
