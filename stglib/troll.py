import warnings
import pandas as pd
import xarray as xr
import numpy as np
import scipy.signal
from .core import utils


def read_aquatroll(filnam, skiprows=69, encoding='utf-8', skipfooter=0):
    return pd.read_csv(filnam,
                       skiprows=skiprows,
                       skipfooter=skipfooter,
                       infer_datetime_format=True,
                       parse_dates=['Date and Time'],
                       encoding=encoding)


def read_aquatroll_header(filnam, encoding='utf-8'):
    with open(filnam, encoding=encoding) as f:
        for line in f.readlines():
            if 'Time Zone:' in line:
                # remove commas and only return the value,
                # not the 'Time Zone: ' part
                return line.replace(',','').strip()[11:]
