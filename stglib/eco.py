#!/usr/bin/env python

from __future__ import division, print_function
import pandas as pd
import xarray as xr

def read_ntu(filnam, spb=False, skiprows=None, skipfooter=0):
    ntu = pd.read_csv(filnam,
                      sep='\t',
                      names=['date', 'time', 'a', 'counts', 'b'],
                      parse_dates=[['date', 'time']],
                      infer_datetime_format=True,
                      engine='python',
                      skiprows=skiprows,
                      skipfooter=skipfooter)

    if spb:
        times = ntu['date_time'].values.reshape((-1, spb))[:, int(spb/2)] # get middle time
        counts = ntu['counts'].values.reshape((-1, spb))
        a = ntu['a'].values.reshape((-1, spb))
        b = ntu['b'].values.reshape((-1, spb))
        sample = range(spb)

        ds = xr.Dataset({'time': ('time', times),
                         'counts': (['time', 'sample'], counts),
                         'a': (['time', 'sample'], a),
                         'b': (['time', 'sample'], b),
                         'sample': ('sample', sample)})
    else:
        times = ntu['date_time']
        counts = ntu['counts']

        ds = xr.Dataset({'time': ('time', times),
                         'counts': ('time', counts)})

    return ds

# TODO: implement read_par()

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Convert raw NTU tab-delimited text files to netCDF')
    parser.add_argument('ntuname', help='raw NTU filename')
    parser.add_argument('--spb', help='samples per burst, default 10', default=10)

    args = parser.parse_args()

    ntu = read_ntu(args.ntuname)

    ntu = pd_to_xr(ntu, spb=args.spb)

    print(ntu)

if __name__ == '__main__':
    main()
