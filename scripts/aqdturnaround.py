#!/usr/bin/env python

import stglib
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np

args = stglib.cmd.aqdturnaround_parser().parse_args()

meta = stglib.aqd.qaqc.read_aqd_hdr(args.basefile)

ds = stglib.aqd.hdr2cdf.load_sen(args.basefile)

ds = stglib.utils.write_metadata(ds, meta)

ds = stglib.aqd.hdr2cdf.load_amp_vel(ds, args.basefile)

T_orig = meta['AQDTransMatrix'].copy()
T = meta['AQDTransMatrix']
if args.orientation == 'DOWN':
    T[1,:] = -T[1,:]
    T[2,:] = -T[2,:]

u, v, w = stglib.aqd.qaqc.coord_transform(ds['VEL1'].values,
                                          ds['VEL2'].values,
                                          ds['VEL3'].values,
                                          ds['Heading'].values,
                                          ds['Pitch'].values,
                                          ds['Roll'].values,
                                          T,
                                          T_orig,
                                          meta['AQDCoordinateSystem'])

ds['U'] = xr.DataArray(u, dims=('time', 'bindist'))
ds['V'] = xr.DataArray(v, dims=('time', 'bindist'))
ds['W'] = xr.DataArray(w, dims=('time', 'bindist'))

ds.attrs['AQDTransMatrix'] = []

ds.to_netcdf(args.basefile + '_turnaround.nc')

plt.figure(figsize=(8.5, 11))
plt.subplot(4, 1, 1)
ds['Pressure'].plot()
plt.title(args.basefile)

plt.subplot(4, 1, 2)
ds['Heading'].plot()
plt.ylim(0, 360)

plt.subplot(4, 1, 3)
ds['Pitch'].plot()
plt.ylim(-15, 15)

plt.subplot(4, 1, 4)
ds['Roll'].plot()
plt.ylim(-15, 15)

plt.savefig(args.basefile + '_aqd_sensor_ts.png')


plt.figure(figsize=(11, 9))
for sp, var in zip([1, 2, 3], ['U', 'V', 'W']):
    plt.subplot(3, 1, sp)
    ds[var].transpose().plot.pcolormesh(
        vmax=np.max([-ds[var].quantile(0.05), ds[var].quantile(0.95)]))
    ds['Pressure'].plot(c='k')
    plt.title(args.basefile)
plt.savefig(args.basefile + '_aqd_velocity_pcolor.png')


plt.figure(figsize=(11, 9))
for sp, var in zip([1, 2, 3], ['1', '2', '3']):
    plt.subplot(3, 1, sp)
    ds['AMP' + var].transpose().plot.pcolormesh(
        vmax=ds['AMP' + var].quantile(0.95))
    ds['Pressure'].plot(c='k')
    plt.title(args.basefile)
plt.savefig(args.basefile + '_aqd_amplitude_pcolor.png')


plt.figure(figsize=(11, 9))
plt.subplot(4, 1, 1)
ds['Pressure'].plot()

for sp in [2, 3, 4]:
    plt.subplot(4, 1, sp)
    for beam in ['U', 'V', 'W']:
        ds[beam].isel(bindist=sp-2).plot()  # second bin
    plt.legend()
    plt.title((args.basefile + ' bin ' + str(sp-2) +
               ' (' + str(ds['bindist'][sp-2].values) + ' m)'))
plt.savefig(args.basefile + '_aqd_velocity_ts.png')

# plt.figure(figsize=(11,9))
# sp = plt.subplot(1, 1, 1, projection='polar')
# sp.set_theta_zero_location('N')
# sp.set_theta_direction(-1)
