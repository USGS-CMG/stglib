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
                       parse_dates={'time': ['Year', 'Month', 'Day', 'Hour', 'Minute', 'Second']},
                       date_parser=qaqc.date_parser)
    df.set_index('time', inplace=True)
    return df

def read_dat(filnam):
    return read_dat_raw(filnam)

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
