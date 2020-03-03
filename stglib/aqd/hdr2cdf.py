from __future__ import division, print_function
import pandas as pd
import xarray as xr
from ..core import utils
from . import qaqc


def prf_to_cdf(metadata):
    """Load a Aquadopp text files and output to netCDF format"""

    # TODO: clock drift code
    # TODO: logmeta code

    basefile = metadata['basefile']

    if 'prefix' in metadata:
        basefile = metadata['prefix'] + basefile

    utils.check_valid_metadata(metadata)

    # get instrument metadata from the HDR file
    instmeta = qaqc.read_aqd_hdr(basefile)

    metadata['instmeta'] = instmeta

    print("Loading ASCII files")

    # Load sensor data
    ds = load_sen(basefile)

    # write out metadata first, then deal exclusively with xarray attrs
    ds = utils.write_metadata(ds, metadata)

    del metadata
    del instmeta

    # Deal with metadata peculiarities
    ds = qaqc.check_attrs(ds)

    ds = qaqc.check_orientation(ds)

    # Load amplitude and velocity data
    ds = load_amp_vel(ds, basefile)

    # Compute time stamps
    ds = utils.shift_time(ds, ds.attrs['AQDAverageInterval']/2)

    if utils.is_cf(ds):
        pass
    else:
        print('about to create epic times')
        ds = utils.create_epic_times(ds)

    # configure file
    if 'prefix' in ds.attrs:
        cdf_filename = ds.attrs['prefix'] + ds.attrs['filename'] + '-raw.cdf'
    else:
        cdf_filename = ds.attrs['filename'] + '-raw.cdf'

    ds = qaqc.update_attrs(ds)

    # need to drop datetime
    ds = ds.drop('datetime')

    ds.to_netcdf(cdf_filename, unlimited_dims=['time'])

    print('Finished writing data to %s' % cdf_filename)

    return ds


def load_sen(basefile):
    """Load data from .sen file"""

    senfile = basefile + '.sen'

    SEN = pd.read_csv(senfile,
                      header=None,
                      delim_whitespace=True,
                      parse_dates={'datetime': [2, 0, 1, 3, 4, 5]},
                      date_parser=qaqc.date_parser,
                      usecols=[0, 1, 2, 3, 4, 5, 8, 10,
                               11, 12, 13, 14, 15, 16])

    # rename columns from numeric to human-readable
    SEN.rename(columns={10: 'Heading',
                        11: 'Pitch',
                        12: 'Roll',
                        13: 'Pressure',
                        14: 'Temperature',
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
        a = pd.read_csv(afile, header=None, delim_whitespace=True)

        if 'bindist' in RAW:
            coords = [RAW['time'], RAW['bindist']]
        else:
            coords = [RAW['time'], RAW.attrs['AQDCCD']]

        RAW['AMP' + str(n)] = xr.DataArray(a,
                                           dims=('time', 'bindist'),
                                           coords=coords)

        vfile = basefile + '.v' + str(n)
        v = pd.read_csv(vfile, header=None, delim_whitespace=True)
        # convert m/s to cm/s
        RAW['VEL' + str(n)] = xr.DataArray(v * 100,
                                           dims=('time', 'bindist'),
                                           coords=coords)

    return RAW
