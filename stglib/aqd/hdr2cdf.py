from __future__ import division, print_function

import warnings
import pandas as pd
import xarray as xr
import numpy as np
from ..core import utils
from . import qaqc

def prf_to_cdf(metadata):
    """Main load file"""

    # TODO: clock drift code
    # TODO: Move time to center of ensemble??
    # TODO: logmeta code

    basefile = metadata['basefile']

    # get instrument metadata from the HDR file
    instmeta = qaqc.read_aqd_hdr(basefile)

    metadata['instmeta'] = instmeta

    print("Loading ASCII files")

    # Load sensor data
    ds = load_sen(basefile)

    # write out metadata first, then deal exclusively with xarray attrs
    ds = utils.write_metadata(ds, metadata)
    ds = utils.write_metadata(ds, metadata['instmeta'])

    del metadata
    del instmeta

    # Deal with metadata peculiarities
    ds = qaqc.check_attrs(ds)

    ds = qaqc.check_orientation(ds)

    # Load amplitude and velocity data
    ds = load_amp_vel(ds, basefile)

    # Compute time stamps
    ds = utils.shift_time(ds, ds.attrs['AQDAverageInterval']/2)

    ds = utils.create_epic_time(ds)

    # configure file
    cdf_filename = ds.attrs['filename'] + '-raw.cdf'

    ds = qaqc.update_attrs(ds)

    # need to drop datetime
    ds = ds.drop('datetime')

    ds.to_netcdf(cdf_filename, unlimited_dims='time')

    print('Finished writing data to %s' % cdf_filename)

    return ds


def load_sen(basefile):
    """Load data from .sen file"""

    senfile = basefile + '.sen'

    # read csv and parse dates
    # https://stackoverflow.com/questions/27112591/parsing-year-month-day-hour-minute-second-in-python
    def parse(year, month, day, hour, minute, second):
        return year + '-' + month + '-' + day + ' ' + hour + ':' + minute + ':' + second

    SEN = pd.read_csv(senfile,
                      header=None,
                      delim_whitespace=True,
                      parse_dates={'datetime': [2, 0, 1, 3, 4, 5]},
                      date_parser=parse,
                      usecols=[0, 1, 2, 3, 4, 5, 8, 10, 11, 12, 13, 14, 15, 16])

    # rename columns from numeric to human-readable
    SEN.rename(columns={10: 'Heading',
                        11: 'Pitch',
                        12: 'Roll',
                        13:'Pressure',
                        14:'Temperature',
                        8: 'Battery'},
               inplace=True)

    # Look for analog data TODO
    SEN.rename(columns={15: 'AnalogInput1'}, inplace=True)
    SEN['AnalogInput1'] = SEN['AnalogInput1'] * 5 / 65535
    SEN.rename(columns={16: 'AnalogInput2'}, inplace=True)
    SEN['AnalogInput2'] = SEN['AnalogInput2'] * 5 / 65535

    # create xarray Dataset
    RAW = xr.Dataset.from_dataframe(SEN)
    RAW = RAW.rename({'index': 'time'})
    RAW['time'] = RAW['datetime']

    return RAW


def load_amp_vel(RAW, basefile):
    """Load amplitude and velocity data from the .aN and .vN files"""

    for n in [1, 2, 3]:
        afile = basefile + '.a' + str(n)
        RAW['AMP' + str(n)] = xr.DataArray(pd.read_csv(afile, header=None, delim_whitespace=True),
            dims=('time', 'bindist'), coords=[RAW['time'], RAW['bindist']])
        # RAW['AMP' + str(n)] = RAW['AMP' + str(n)].rename({'dim_0': 'time'})
        vfile = basefile + '.v' + str(n)
        v = pd.read_csv(vfile, header=None, delim_whitespace=True)
        # convert to cm/s
        RAW['VEL' + str(n)] = xr.DataArray(v * 100, dims=('time', 'bindist'), coords=[RAW['time'], RAW['bindist']])

    return RAW


def insert_fill_values(RAW):
    """Insert fill values for nans"""

    print("Inserting fill values")
    for k in RAW:
        if k not in ['instmeta', 'time', 'time2', 'datetime'] and np.max(np.shape(RAW[k])) == np.max(np.shape(RAW['jd'])):
            nanind = np.where(np.isnan(RAW[k]))
            RAW[k][nanind] = 1e35

    return RAW
