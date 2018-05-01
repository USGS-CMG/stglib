from __future__ import division, print_function
import scipy.signal as spsig
import numpy as np


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
    return 1/Pnn['frequency'][Pnn.fillna(0).argmax(dim='frequency')].values


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
