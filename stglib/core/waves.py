from __future__ import division, print_function
import scipy.signal as spsig
import numpy as np
import xarray as xr

def make_waves_ds(ds, noise=0.9):

    print('Computing waves statistics')

    f, Pxx = pressure_spectra(ds['P_1ac'],
                                    fs=1/ds.attrs['sample_interval'])

    z = ds.attrs['initial_instrument_height']
    h = ds['P_1ac'].mean(dim='sample') + z

    k = np.asarray(
        [qkfs(2*np.pi/(1/f), x) for x in h.values])

    Kp = transfer_function(k, h, z)
    Pnn = elevation_spectra(Pxx, Kp)

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
        *[define_cutoff(f, x, noise=0.75) for x in spec['Pxx'].values])
    spec['tailind'] = xr.DataArray(np.asarray(tailind), dims='time')
    spec['noisecutind'] = xr.DataArray(np.asarray(noisecutind), dims='time')
    spec['fpeakcutind'] = xr.DataArray(np.asarray(fpeakcutind), dims='time')
    thetail = [make_tail(
        spec['frequency'],
        spec['Pnn'][burst, :],
        spec['tailind'][burst].values)
        for burst in range(len(spec['time']))]
    spec['pspec'] = xr.DataArray(thetail, dims=('time', 'frequency'))
    spec['m0'] = xr.DataArray(
        make_m0(spec['frequency'], spec['pspec']),
        dims='time')
    spec['m2'] = xr.DataArray(
        make_m2(spec['frequency'], spec['pspec']),
        dims='time')
    spec['wh_4061'] = xr.DataArray(
        make_Hs(spec['m0']), dims='time')
    spec['wp_4060'] = xr.DataArray(
        make_Tm(spec['m0'], spec['m2']), dims='time')
    spec['wp_peak'] = xr.DataArray(make_Tp(spec['pspec']), dims='time')
    spec['kh'] = xr.DataArray(k, dims=('time', 'frequency'))

    return spec


def pressure_spectra(x, fs=1.0, window='hanning', nperseg=256, **kwargs):
    """Compute pressure spectral density using Welch's method

    Parameters
    ----------
    x : array_like
        Time-series of pressure data
    fs : float, optional
        Sampling frequency (Hz)
    window : str, optional
        Window, default 'hanning'
    nperseg : int, optional
        Length of each segment, default 256
    **kwargs
        Arbitrary keyword arguments passed to scipy.signal.welch

    Returns
    -------
    f : ndarray
        Array of sample frequencies
    Pxx : ndarray
        Power spectral density of pressure data
    """
    f, Pxx = spsig.welch(x, fs=fs, window=window, nperseg=nperseg, **kwargs)
    return f, Pxx


def elevation_spectra(Pxx, Kp):
    """Compute elevation spectra using linear wave theory and transfer function
    """
    return Pxx / (Kp**2)


def transfer_function(k, h, z):
    """Compute pressure transfer function

    Parameters
    ----------
    k : float
        Wavenumber
    h : float
        Water depth [m]
    z: float
        Height of pressure sensor above bottom [m]

    Returns
    -------
    Kp : float
        Presssure transfer function
    """
    Kp = np.cosh(k*z)/np.cosh(k*np.expand_dims(h, 1))

    # set Kp nans at 0 frequency to 1
    Kp[np.isnan(k)] = 1

    return Kp


def define_cutoff(f, Pxx, noise=0.9):
    noisecut = 12*np.mean(Pxx[f >= noise*f[-1]])
    tmp = np.where(Pxx <= noisecut)[0]
    # sometimes the first value is less than the noise floor (not sure why)
    if 0 in tmp:  # it has chosen the first value, which we want to ignore
        noisecutind = tmp[1] - 1
    else:
        noisecutind = tmp[0] - 1
    fpeakcut = 1.1*f[np.argmax(Pxx)]
    fpeakcutind = np.searchsorted(f, fpeakcut)

    if noisecutind > fpeakcutind:
        tailind = noisecutind
    else:
        tailind = np.nan
    return tailind, noisecut, noisecutind, fpeakcutind


def make_tail(f, Pnn, tailind):
    ti = tailind.astype(int)
    tail = np.ones_like(f)
    if np.isnan(tailind):
        return np.ones_like(f) * np.nan
    else:
        tail[ti:] = Pnn[ti] * (f[ti:]/f[ti])**-4
        return np.hstack((Pnn[:ti], tail[ti:]))


def make_m0(f, Pnn):
    return np.trapz(Pnn, x=f)


def make_m2(f, Pnn):
    return np.trapz(Pnn * f**2, x=f)


def make_Hs(m0):
    return 4*np.sqrt(m0)


def make_Tm(m0, m2):
    return np.sqrt(m0/m2)


def make_Tp(Pnn):
    # ensure we don't return 0 frequency as a peak period
    fp = Pnn['frequency'][Pnn.fillna(0).argmax(dim='frequency')].values
    fp[fp == 0] = np.nan
    return 1/fp


def qkfs(omega, h):
    """
    Modified from Wiberg & Sherwood 2009; only does 3 iterations.
    Returns only k, not kh

    k = qkfs(omega, h)
    """

    g = 9.81
    x = omega**2 * h / g
    y = np.sqrt(x) * (x < 1) + x * (x >= 1)

    t = np.tanh(y)
    y = y - ((y*t-x) / (t + y * (1-t**2)))
    t = np.tanh(y)
    y = y - ((y*t-x) / (t + y * (1-t**2)))
    t = np.tanh(y)
    y = y - ((y*t-x) / (t + y * (1-t**2)))

    return y/h