from __future__ import division, print_function

import pandas as pd
import xarray as xr
import numpy as np
from ..core import utils
from . import qaqc


def wad_to_cdf(metadata, writefile=True):
    """Load Aquadopp waves data and create raw netCDF file

    Parameters
    ----------
    metadata : dict
        Dictionary of required metadata
    writefile : bool, optional
        Flag to write raw .cdf file. Default True

    Returns
    -------
    xarray.Dataset
        Raw waves data in an xarray Dataset

    """

    basefile = metadata['basefile']

    # get instrument metadata from the HDR file
    instmeta = qaqc.read_aqd_hdr(basefile)

    metadata['instmeta'] = instmeta

    ds = load_whd(metadata)

    # write out metadata first, then deal exclusively with xarray attrs
    ds = utils.write_metadata(ds, metadata)
    ds = utils.write_metadata(ds, metadata['instmeta'])

    del metadata
    del instmeta

    ds = load_wad(ds)

    # Deal with metadata peculiarities
    ds = qaqc.check_attrs(ds, waves=True)

    ds.attrs['center_first_bin'] = ds['cellpos'][0].values

    print('BIN SIZE:', ds.attrs['bin_size'])

    ds = qaqc.check_orientation(ds, waves=True)

    # Compute time stamps
    fs = float(ds.attrs['WaveSampleRate'].split()[0])
    ds.attrs['sample_interval'] = 1/fs
    ds.attrs['samples_per_burst'] = ds.attrs['WaveNumberOfSamples']

    if utils.is_cf(ds):
        pass
    else:
        print('about to create epic times')
        ds = utils.create_epic_times(ds, waves=True)

        ds = utils.create_2d_time(ds)

    ds = utils.ds_coord_no_fillvalue(ds)

    ds = qaqc.update_attrs(ds, waves=True)

    # need to drop datetime
    ds = ds.drop('datetime')

    if writefile:
        cdf_filename = ds.attrs['filename'] + 'wvs-raw.cdf'
        ds.to_netcdf(cdf_filename)
        print('Finished writing data to %s' % cdf_filename)

    return ds


def load_whd(metadata):
    """Load data from .whd file"""

    whdfile = metadata['basefile'] + '.whd'

    WHD = pd.read_csv(
        whdfile,
        header=None,
        delim_whitespace=True,
        parse_dates={'datetime': [2, 0, 1, 3, 4, 5]},
        date_parser=qaqc.date_parser,
        usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                 10, 11, 12, 13, 14, 16, 17, 18, 19, 20])

    # rename columns from numeric to human-readable
    WHD.rename(columns={6: 'burst',
                        7: 'nrecs',
                        8: 'cellpos',
                        9: 'Battery',
                        10: 'soundspeed',
                        11: 'Heading',
                        12: 'Pitch',
                        13: 'Roll',
                        14: 'minpressure',
                        16: 'Temperature',
                        17: 'cellsize',
                        18: 'avgamp1',
                        19: 'avgamp2',
                        20: 'avgamp3'},
               inplace=True)

    ds = xr.Dataset.from_dataframe(WHD)
    ds = ds.rename({'index': 'time'})
    ds['time'] = ds['datetime']
    ds = ds.drop(['minpressure', 'cellsize', 'nrecs', 'soundspeed'])

    return ds


def load_wad(ds):

    wadfile = ds.attrs['basefile'] + '.wad'
    print('Loading wave data from ' + wadfile + '; this may take some time')
    # pd.read_csv is ~10x faster than np.loadtxt or np.genfromtxt
    WAD = pd.read_csv(wadfile, header=None, delim_whitespace=True).values

    r, c = np.shape(WAD)
    print(wadfile + ' has ' + str(r) + ' rows and ' + str(c) + ' columns')
    if 'num_wave_bursts' in ds.attrs: # we can override the number of samples if need be
        print('Overriding number of samples using attr num_wave_bursts of {}'.format(ds.attrs['num_wave_bursts']))
        ds = ds.sel(time=ds.time[0:ds.attrs['num_wave_bursts']])
        nburst = ds.attrs['num_wave_bursts']
    elif r % ds.attrs['WaveNumberOfSamples']:
        print('Number of rows read is not a multiple of %d. Truncating data '
              'to last full burst' %
              ds.attrs['WaveNumberOfSamples'])
        ds = ds.sel(time=ds.time[0:int(np.floor(r / ds.attrs['WaveNumberOfSamples']))])
        nburst = int(np.floor(r/ds.attrs['WaveNumberOfSamples']))
    else:
        nburst = int(np.floor(r/ds.attrs['WaveNumberOfSamples']))
    nsamps = int(nburst * ds.attrs['WaveNumberOfSamples'])
    wavensamps = int(ds.attrs['WaveNumberOfSamples'])
    print('Metadata reports ' + str(nburst) + ' bursts, ' + str(nsamps) +
          ' samples, ' + str(wavensamps) + ' samples per burst')

    thevars = ['Pressure', 'VEL1', 'VEL2', 'VEL3', 'AMP1', 'AMP2', 'AMP3']
    thecols = [2, 5, 6, 7, 9, 10, 11]
    for var, n in zip(thevars, thecols):
        ds[var] = xr.DataArray(
            np.reshape(WAD[0:nsamps, n], (nburst, wavensamps)),
            dims=('time', 'sample'))

    # convert m/s to mm/s
    for var in ['VEL1', 'VEL2', 'VEL3']:
        ds[var] = ds[var] * 1000

    print('Done loading ' + wadfile)

    return ds
