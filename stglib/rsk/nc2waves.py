from __future__ import division, print_function
import xarray as xr
import numpy as np
from ..core import utils, waves


def nc_to_waves(nc_filename):

    ds = xr.open_dataset(nc_filename, autoclose=True, decode_times=False)

    ds = utils.epic_to_cf_time(ds)

    ds = utils.create_epic_time(ds)

    spec = make_waves_ds(ds)

    for k in ['wp_peak', 'wh_4061', 'wp_4060', 'pspec']:
        ds[k] = spec[k]

    ds = utils.create_water_depth(ds)

    ds = ds.drop(['P_1', 'P_1ac', 'sample'])

    ds = utils.trim_max_wp(ds)

    ds = utils.trim_min_wh(ds)

    ds = utils.trim_max_wh(ds)

    ds = utils.trim_wp_ratio(ds)

    # Add attrs
    ds = utils.ds_add_attrs(ds)

    # Reshape and associate dimensions with lat/lon
    for var in ['wp_peak', 'wh_4061', 'wp_4060', 'pspec']:
        if var in ds:
            ds = utils.add_lat_lon(ds, var)

    # assign min/max (need to do this after trimming):
    ds = utils.add_min_max(ds)

    ds = utils.ds_add_diwasp_history(ds)

    nc_filename = ds.attrs['filename'] + 's-a.nc'

    ds = utils.rename_time(ds)

    ds.to_netcdf(nc_filename)

    return ds


def make_waves_ds(ds, noise=0.9):

    print('Computing waves statistics')

    f, Pxx = waves.pressure_spectra(ds['P_1ac'],
                                    fs=1/ds.attrs['sample_interval'])

    z = ds.attrs['initial_instrument_height']
    h = ds['P_1ac'].mean(dim='sample') + z

    k = np.asarray(
        [waves.qkfs(2*np.pi/(1/f), x) for x in h.values])

    Kp = waves.transfer_function(k, h, z)
    Pnn = waves.elevation_spectra(Pxx, Kp)

    spec = xr.Dataset()
    # spec['time'] = xr.DataArray(dw1076['b']['P_1ac'].time, dims='time')
    # spec['frequency'] = xr.DataArray(f, dims='frequency')
    spec['Pnn'] = xr.DataArray(Pnn,
                               dims=('time', 'frequency'),
                               coords=(ds['time'], f))
    spec['Pxx'] = xr.DataArray(Pxx,
                               dims=('time', 'frequency'),
                               coords=(ds['time'], f))
    tailind, noisecut, noisecutind, fpeakcutind = zip(
        *[waves.define_cutoff(f, x, noise=0.75) for x in spec['Pxx'].values])
    spec['tailind'] = xr.DataArray(np.asarray(tailind), dims='time')
    spec['noisecutind'] = xr.DataArray(np.asarray(noisecutind), dims='time')
    spec['fpeakcutind'] = xr.DataArray(np.asarray(fpeakcutind), dims='time')
    thetail = [waves.make_tail(
        spec['frequency'],
        spec['Pnn'][burst, :],
        spec['tailind'][burst].values)
        for burst in range(len(spec['time']))]
    spec['pspec'] = xr.DataArray(thetail, dims=('time', 'frequency'))
    spec['m0'] = xr.DataArray(
        waves.make_m0(spec['frequency'], spec['pspec']),
        dims='time')
    spec['m2'] = xr.DataArray(
        waves.make_m2(spec['frequency'], spec['pspec']),
        dims='time')
    spec['wh_4061'] = xr.DataArray(
        waves.make_Hs(spec['m0']), dims='time')
    spec['wp_4060'] = xr.DataArray(
        waves.make_Tm(spec['m0'], spec['m2']), dims='time')
    spec['wp_peak'] = xr.DataArray(waves.make_Tp(spec['pspec']), dims='time')
    spec['kh'] = xr.DataArray(k, dims=('time', 'frequency'))

    return spec
