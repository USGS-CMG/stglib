#!/usr/bin/env python

from __future__ import division, print_function
import sys
sys.path.insert(0, '/Users/dnowacki/Documents/aqdlib')
from aqdlib.aqdhdr2cdf import compute_time, read_aqd_hdr, check_orientation, check_metadata, write_metadata, update_attrs
import numpy as np
import pandas as pd
import xarray as xr

def wad_to_cdf(metadata):
    """Main waves load file"""

    basefile = metadata['basefile']

    # get instrument metadata from the HDR file
    instmeta = read_aqd_hdr(basefile)

    metadata['instmeta'] = instmeta

    RAW = load_whd(metadata)

    RAW = load_wad(RAW, metadata)

    # Deal with metadata peculiarities
    metadata = check_metadata(metadata, waves=True)

    metadata['center_first_bin'] = RAW['cellpos'][0].values

    print('BIN SIZE:', metadata['bin_size'])

    RAW = check_orientation(RAW, metadata, waves=True)

    # Compute time stamps
    RAW = compute_time(RAW, metadata, waves=True)

    # configure file
    cdf_filename = metadata['filename'] + 'wvs-raw.cdf' # TODO: fix the path

    # write out metadata
    RAW = write_metadata(RAW, metadata)
    RAW = write_metadata(RAW, metadata['instmeta'])

    update_attrs(cdf_filename, RAW, metadata, waves=True)

    # need to drop datetime
    RAW = RAW.drop('datetime')

    RAW.to_netcdf(cdf_filename, unlimited_dims='time')

    print('Finished writing data to %s' % cdf_filename)

    return RAW, metadata

def load_whd(metadata):
    """Load data from .whd file"""

    whdfile = metadata['basefile'] + '.whd'

    # read csv and parse dates
    # https://stackoverflow.com/questions/27112591/parsing-year-month-day-hour-minute-second-in-python
    def parse(year, month, day, hour, minute, second):
        return year + '-' + month + '-' + day + ' ' + hour + ':' + minute + ':' + second

    WHD = pd.read_csv(whdfile,
        header=None,
        delim_whitespace=True,
        parse_dates={'datetime': [2, 0, 1, 3, 4, 5]},
        date_parser=parse,
        usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20])

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

    RAW = xr.Dataset.from_dataframe(WHD)
    RAW = RAW.rename({'index': 'time'})
    RAW['time'] = RAW['datetime']

    return RAW

def load_wad(RAW, metadata):

    wadfile = metadata['basefile'] + '.wad'
    print('Loading wave data from ' + wadfile + '; this may take some time')
    # pd.read_csv is ~10x faster than np.loadtxt or np.genfromtxt
    WAD = pd.read_csv(wadfile, header=None, delim_whitespace=True).values

    r, c = np.shape(WAD)
    print(r, c)
    nburst = int(np.floor(r/metadata['instmeta']['WaveNumberOfSamples']))
    nsamps = int(nburst * metadata['instmeta']['WaveNumberOfSamples'])
    wavensamps = int(metadata['instmeta']['WaveNumberOfSamples'])
    print(nburst, nsamps, wavensamps)

    samples = np.arange(wavensamps)

    RAW['sample'] = xr.DataArray(samples, dims=('sample'), name='sample')

    RAW['Pressure'] = xr.DataArray(np.reshape(WAD[0:nsamps, 2], (nburst, wavensamps)), dims=('time', 'sample'))
    RAW['VEL1'] = xr.DataArray(np.reshape(WAD[0:nsamps, 5], (nburst, wavensamps)), dims=('time', 'sample'))
    RAW['VEL2'] = xr.DataArray(np.reshape(WAD[0:nsamps, 6], (nburst, wavensamps)), dims=('time', 'sample'))
    RAW['VEL3'] = xr.DataArray(np.reshape(WAD[0:nsamps, 7], (nburst, wavensamps)), dims=('time', 'sample'))
    RAW['AMP1'] = xr.DataArray(np.reshape(WAD[0:nsamps, 9], (nburst, wavensamps)), dims=('time', 'sample'))
    RAW['AMP2'] = xr.DataArray(np.reshape(WAD[0:nsamps, 10], (nburst, wavensamps)), dims=('time', 'sample'))
    RAW['AMP3'] = xr.DataArray(np.reshape(WAD[0:nsamps, 11], (nburst, wavensamps)), dims=('time', 'sample'))

    # convert to cm/s
    for n in [1, 2, 3]:
        RAW['VEL' + str(n)] = RAW['VEL' + str(n)] * 100

    print('Done loading ' + wadfile)

    return RAW

def main():
    import sys
    sys.path.insert(0, '/Users/dnowacki/Documents/aqdlib')
    import aqdlib
    import argparse
    import yaml

    parser = argparse.ArgumentParser(description='Convert Aquadopp .wad wave files to raw .cdf format. Run this script from the directory containing Aquadopp files')
    # parser.add_argument('basename', help='base name (without extension) of the Aquadopp text files')
    parser.add_argument('gatts', help='path to global attributes file (gatts formatted)')
    parser.add_argument('config', help='path to ancillary config file (YAML formatted)')
    # parser.add_argument('-v', '--verbose', help='increase output verbosity', action='store_true')

    args = parser.parse_args()

    # initialize metadata from the globalatts file
    metadata = aqdlib.read_globalatts(args.gatts)

    # Add additional metadata from metadata config file
    config = yaml.safe_load(open(args.config))

    for k in config:
        metadata[k] = config[k]

    wad_to_cdf(metadata)

if __name__ == '__main__':
    main()
