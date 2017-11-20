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

# # This is when it was sampling properly
# ds = read_ntu(filnam, spb=10, skiprows=151, skipfooter=107079)
# # This is when it was sampling continuously
# ds2 = read_ntu(filnam, skiprows=1161, skipfooter=1)
#
# print(ds)
# print(ds2)
# # %%
# plt.figure(figsize=(11,8.5))
# plt.plot(ds2['time'], ds2['counts'])
# plt.show()
#
# # %%
#
# plt.figure(figsize=(11,8.5))
# plt.plot(ds['time'], ds['counts'].mean(dim='sample'))
# plt.ylabel('NTU')
# # plt.savefig('/Volumes/Backstaff/field/nrpp/gate/1100ntu.pdf')
# plt.show()

# %%

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
