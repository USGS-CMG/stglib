from __future__ import division, print_function
import scipy.signal as spsig
import numpy as np
import xarray as xr


def make_waves_ds(ds, noise=0.75):

    print('Computing waves statistics')

    f, Pxx = pressure_spectra(ds['P_1ac'].squeeze(),
                              fs=1/ds.attrs['sample_interval'])

    z = ds.attrs['initial_instrument_height']
    h = ds['P_1ac'].squeeze().mean(dim='sample') + z

    k = np.asarray(
        [qkfs(2*np.pi*f, x) for x in h.values])

    Kp = transfer_function(k, h, z)
    Pnn = elevation_spectra(Pxx, Kp)

    spec = xr.Dataset()

    spec['Pnn'] = xr.DataArray(Pnn,
                               dims=('time', 'frequency'),
                               coords=(ds['time'], f))
    spec['Pxx'] = xr.DataArray(Pxx,
                               dims=('time', 'frequency'),
                               coords=(ds['time'], f))
    spec['Kp'] = xr.DataArray(Kp,
                              dims=('time', 'frequency'),
                              coords=(ds['time'], f))
    tailind, noisecutind, fpeakcutind, Kpcutind = zip(
        *[define_cutoff(f, x, y, noise=noise) for x, y in
            zip(spec['Pxx'].values, spec['Kp'].values)])
    spec['tailind'] = xr.DataArray(np.asarray(tailind), dims='time')
    spec['noisecutind'] = xr.DataArray(np.asarray(noisecutind), dims='time')
    spec['fpeakcutind'] = xr.DataArray(np.asarray(fpeakcutind), dims='time')
    spec['Kpcutind'] = xr.DataArray(np.asarray(Kpcutind), dims='time')
    thetail = [make_tail(
        spec['frequency'],
        spec['Pnn'][burst, :],
        spec['tailind'][burst].values)
        for burst in range(len(spec['time']))]
    spec['pspec'] = xr.DataArray(thetail, dims=('time', 'frequency'))
    spec['m0'] = xr.DataArray(
        make_moment(spec['frequency'], spec['pspec'], 0),
        dims='time')
    spec['m2'] = xr.DataArray(
        make_moment(spec['frequency'], spec['pspec'], 2),
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


def define_cutoff(f, Pxx, Kp, noise=0.9):
    """Define cutoff based on Jones & Monismith (2007) values for the
    noise floor cutoff (12*noise), f_peak cutoff (1.1*f_peak), along with
    a cutoff based on 0.1 K_p (pressure transfer function)

    Parameters
    ----------
    f : array_like
        Frequencies
    Pxx : array_like
        Untransformed pressure spectra
    Kp : array_like
        Pressure transfer function
    noise : float, optional
        Fractional frequency above which we consider data to be noise.
        Default 0.9 (i.e. top 10% of frequencies considered noise)

    Returns
    -------
    tailind : int
        Index for where to start applying the f^-4 tail
    noisecutind : int
        Index for noise cutoff
    fpeakcutind : int
        Index for f_peak cutoff
    Kpcutind : int
        Index for the K_p cutoff

    References
    ----------
    Jones, N. L., & Monismith, S. G. (2007). Measuring short-period wind waves
    in a tidally forced environment with a subsurface pressure gauge.
    Limnology and Oceanography: Methods, 5, 317–327.
    http://doi.org/10.4319/lom.2007.5.317
    """

    noisecut = 12*np.mean(Pxx[f >= noise*f[-1]])
    # Look for above cut-off freqs and take highest
    tmp = np.where(Pxx > noisecut)[0]
    if len(tmp) == 0:
        noisecutind = 0
    else:
        noisecutind = tmp[-1]

    fpeakcut = 1.1*f[np.argmax(Pxx)]
    fpeakcutind = np.searchsorted(f, fpeakcut)  # cutoff based on 1.1*fp
    Kpcutind = np.argmax(Kp <= 0.1)  # cutoff based on Kp<=0.1

    if (noisecutind > fpeakcutind) and (noisecutind <= Kpcutind):
        tailind = noisecutind
    elif (noisecutind > fpeakcutind) and (noisecutind > Kpcutind):
        tailind = Kpcutind
    else:
        tailind = np.nan
    return tailind, noisecutind, fpeakcutind, Kpcutind


def make_tail(f, Pnn, tailind):
    """Make f^-4 tail following Jones & Monismith (2007)

    Parameters
    ----------
    f : array_like
        Frequencies
    Pnn : array_like
        Spectra
    tailind : int
        Index for where to start applying the f^-4 tail

    Returns
    -------
    array_like
        Spectra with f^-4 tail applied above tailind

    References
    ----------
    Jones, N. L., & Monismith, S. G. (2007). Measuring short-period wind waves
    in a tidally forced environment with a subsurface pressure gauge.
    Limnology and Oceanography: Methods, 5, 317–327.
    http://doi.org/10.4319/lom.2007.5.317
    """
    ti = tailind.astype(int)
    tail = np.ones_like(f)
    if np.isnan(tailind):
        return np.ones_like(f) * np.nan
    else:
        tail[ti:] = Pnn[ti] * (f[ti:]/f[ti])**-4
        return np.hstack((Pnn[:ti], tail[ti:]))


def make_mwd(freqs, dirs, dspec):
    """Create mean wave direction (EPIC 4062) variable"""

    Sxsin = dspec * np.expand_dims(np.sin(np.deg2rad(dirs)), axis=1)
    Sxcos = dspec * np.expand_dims(np.cos(np.deg2rad(dirs)), axis=1)

    Dnum = np.trapz(np.trapz(Sxsin, x=freqs), x=dirs)
    Ddnom = np.trapz(np.trapz(Sxcos, x=freqs), x=dirs)

    Dm = np.rad2deg(np.arctan(np.abs(Dnum/Ddnom)))

    Dm[(Dnum > 0) & (Ddnom < 0)] = 180 - Dm[(Dnum > 0) & (Ddnom < 0)]
    Dm[(Dnum < 0) & (Ddnom < 0)] = 180 + Dm[(Dnum < 0) & (Ddnom < 0)]
    Dm[(Dnum < 0) & (Ddnom > 0)] = 360 - Dm[(Dnum < 0) & (Ddnom > 0)]

    return Dm


def make_moment(f, Pnn, n):
    """Compute nth moment (m0, m1, m2, etc.) of power spectra"""
    return np.trapz(Pnn * f**n, x=f)


def make_Hs(m0):
    return 4*np.sqrt(m0)


def make_Tm(m0, m2):
    return np.sqrt(m0/m2)


def make_Tp(Pnn):
    # ensure we don't return 0 frequency as a peak period
    fp = Pnn['frequency'][Pnn.fillna(0).argmax(dim='frequency')].values
    fp[fp == 0] = np.nan
    return 1/fp


def polar2compass(polar):
    """Convert polar directions (starts at positive x, increases
    counterclockwise) to compass directions (starts at positive y, increases
    clockwise)

    http://nautilus.baruch.sc.edu/CSV/explain_coord_sys_diagram.pdf
    """
    comp = -np.atleast_1d(polar) + 90
    comp[comp < 0] = comp[comp < 0] + 360

    return comp


def to2from(todir):
    """Convert "bearing to" directions to "bearing from" directions
    (helpful for waves from/to conventions)
    """

    fromdir = np.atleast_1d(todir) - 180
    fromdir[fromdir < 0] = fromdir[fromdir < 0] + 360

    return fromdir


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
