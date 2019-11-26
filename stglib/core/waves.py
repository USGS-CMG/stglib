from __future__ import division, print_function
import scipy.signal as spsig
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

def make_waves_ds(ds, noise=0.75):

    print('Computing waves statistics')
    if 'P_1ac' in ds:
        presvar = 'P_1ac'
    else:
        presvar = 'P_1'
        
    f, Pxx = pressure_spectra(ds[presvar].squeeze(),
                              fs=1/ds.attrs['sample_interval'])

    z = ds.attrs['initial_instrument_height']
    h = ds[presvar].squeeze().mean(dim='sample') + z

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


def puv_quick(pressure, u, v, depth, height_of_pressure, height_of_velocity, sampling_frequency, fft_length=512,
              rho=1025., first_frequency_cutoff=1 / 50, infra_gravity_cutoff=0.05, last_frequency_cutoff=1 / 5,
              fft_window_type='hanning', show_diagnostic_plot=False, check_variances=False, variance_error=0.0,
              overlap_length='default'):
    """
    Determine wave heights from pressure, east_velocity, v velocity data

    Parameters
    ----------
    pressure : array_like
        pressure (dbar)
    u :  array_like
        u velocities (m/s)
    v :  array_like
        v velocities (m/s)
    depth : float
        Average water depth (m, positive number)
    height_of_pressure : float
        Height of pressure sensor off bottom (m, positive number)
    height_of_velocity : float
        Height of velocity sensor off bottom (m, positive number)
    sampling_frequency : float
        Hz
    fft_length : int
        Length of data to window and process
    rho : float
        Water density (kg/m^3)
    fft_window_type : str
        Data fft_window for spectral calculation, per scipy signal package
    first_frequency_cutoff : float
        Low-frequency cutoff for wave motions
    infra_gravity_cutoff : float
        Infra-gravity wave frequency cutoff
    last_frequency_cutoff : float
        High-frequency cutoff for wave motions
    show_diagnostic_plot : bool
        print plots and other checks
    check_variances : bool
        test to see if variance is preserved in calculations
    variance_error : float
        tolerance for variance preservation, in percent
    overlap_length : str "default" or int length, default will result in fft_length / 2

    Returns
    -------
    dict::
        'Hrmsp': Hrms (=Hmo) from pressure
        'Hrmsu': Hrms from u,v
        'ubr': Representative orbital velocity amplitude in freq. band
            ( first_frequency_cutoff <= f <= last_frequency_cutoff ) (m/s)
        'omegar': Representative orbital velocity (radian frequency)
        'Tr': Representative orbital velocity period (s)
        'Tpp': Peak period from pressure (s)
        'Tpu': Peak period from velocity (s)
        'phir': Representative orbital velocity direction (angles from x-axis, positive ccw)
        'azr': Representative orb. velocity direction (deg; geographic azimuth; ambiguous =/- 180 degrees)
        'ublo': ubr in freq. band (f <= first_frequency_cutoff) (m/s)
        'ubhi': ubr in freq. band (f >= last_frequency_cutoff) (m/s)
        'ubig': ubr in infra-gravity freq. band (first_frequency_cutoff f <= 1/20) (m/s)
        'figure': figure handle
        'axis': axis handle
        'variance_test_passed': True if passing

    References
    ----------
    Madsen (1994) Coastal Engineering 1994, Proc., 24th, Intl. Conf., Coastal Eng. Res. Council / ASCE. pp.384-398.
        (esp. p. 395)
    Thorton & Guza

    Acknowledgements
    ----------------
    converted to python and updated by Marinna Martini from Chris Sherwood's puvq.m.
    puvq.m also had contributions from Laura Landerman and Patrick Dickudt
    """

    gravity = 9.81  # m/s^2
    if fft_window_type is 'hanning':
        fft_window_type = 'hann'  # this is just the way scipy signal likes it
    if overlap_length is "default":
        overlap_length = int(np.floor(fft_length / 2))

    pressure = spsig.detrend(pressure)
    u = spsig.detrend(u)
    v = spsig.detrend(v)

    # compute wave height from velocities

    # Determine velocity spectra for u and v
    [frequencies, Gpp] = spsig.welch(rho * gravity * pressure, fs=sampling_frequency, window=fft_window_type,
                                     nperseg=fft_length, noverlap=overlap_length)

    df = frequencies[2] - frequencies[1]
    [_, Guu] = spsig.welch(u, fs=sampling_frequency, window=fft_window_type, nperseg=fft_length,
                           noverlap=overlap_length)
    [_, Gvv] = spsig.welch(v, fs=sampling_frequency, window=fft_window_type, nperseg=fft_length,
                           noverlap=overlap_length)

    # determine wave number
    omega = np.array([2 * np.pi * x for x in frequencies])  # omega must be numpy array for qkfs
    # catch numpy errors
    np.seterr(divide='ignore', invalid='ignore')
    k = qkfs(omega, float(depth))  # make sure it is float, or qkfs will bomb
    np.seterr(divide=None, invalid=None)

    # compute linear wave transfer function
    kh = k * depth
    kzp = k * height_of_pressure
    kzuv = k * height_of_velocity
    nf = len(omega)
    Hp = np.ones(nf)
    Huv = np.ones(nf)

    # change wavenumber at 0 Hz to 1 to avoid divide by zero
    i = np.array(range(nf))  # this is an index, thus needs to start at first element, in this case 0
    # for some reason in the MATLAB version CRS tests omega for nans instead of k.
    # Here we test k also because that's where the nans show up
    if np.isnan(omega[0]) or np.isnan(k[0]) or (omega[0] <= 0):  # 0 Hz is the first element
        i = i[1:]
        Hp[0] = 1
        Huv[0] = 1

    Hp[i] = rho * gravity * (np.cosh(kzp[i]) / np.cosh(kh[i]))
    Huv[i] = omega[i] * (np.cosh(kzuv[i]) / np.sinh(kh[i]))

    # combine horizontal velocity spectra
    Guv = Guu + Gvv

    # create cut off frequency, so noise is not magnified
    # at least in first testing, subtracting 1 here got closer to the intended freq. cutoff value
    ff = np.argmax(frequencies > first_frequency_cutoff) - 1
    lf = np.argmax(frequencies > last_frequency_cutoff)

    # Determine wave height for velocity spectra
    Snp = Gpp[ff:lf] / (Hp[ff:lf] ** 2)
    Snu = Guv[ff:lf] / (Huv[ff:lf] ** 2)
    fclip = frequencies[ff:lf]

    # Determine rms wave height (multiply by another sqrt(2) for Hs)
    # Thornton and Guza say Hrms = sqrt(8 mo)
    Hrmsu = 2 * np.sqrt(2 * np.sum(Snu * df))
    Hrmsp = 2 * np.sqrt(2 * np.sum(Snp * df))

    # These are representative orbital velocities for w-c calculations,
    # according to Madsen (1994) Coastal Engineering 1994, Proc., 24th
    # Intl. Conf., Coastal Eng. Res. Council / ASCE. pp.384-398.
    # (esp. p. 395)
    ubr = np.sqrt(2 * np.sum(Guv[ff:lf] * df))
    ubr_check = np.sqrt(2 * np.var(u) + 2 * np.var(v))
    omegar = np.sum(omega[ff:lf] * Guv[ff:lf] * df) / np.sum(Guv[ff:lf] * df)
    Tr = 2 * np.pi / omegar

    if len(np.where(np.isnan(Snp))) > 0 | len(np.where(Snp == 0)) > 0:
        Tpp = np.nan
    else:
        jpeak = np.argmax(Snp)  # index location of the maximum value
        Tpp = 1 / fclip[jpeak]

    if len(np.where(np.isnan(Snu))) > 0 | len(np.where(Snu == 0)) > 0:
        Tpu = np.nan
    else:
        jpeak = np.argmax(Snu)
        Tpu = 1 / fclip[jpeak]

    # phi is angle wrt to x axis; this assumes Guu is in x direction
    # phir = atan2( sum(Guu(ff:lf)*df), sum(Gvv(ff:lf)*df) );

    # this is the line changed on 6/24/03 - I still think it is wrong (CRS)
    # phir = atan2( sum(Gvv(ff:lf)*df), sum(Guu(ff:lf)*df) );

    # This is Jessie's replacement for direction
    # 12/08 Jessie notes that Madsen uses velocity and suggests
    # Suu = sqrt(Guu);
    # Svv = sqrt(Gvv);
    # Suv = sqrt(Guv);
    # but I (CRS) think eqn. 24 is based on u^2, so following is ok:
    rr = np.corrcoef(u, v)
    ortest = np.sign(rr[1][0])
    phir = np.arctan2(ortest * np.sum(Gvv[ff:lf] * df), np.sum(Guu[ff:lf] * df))

    # convert to degrees; convert to geographic azimuth (0-360, 0=north)
    azr = 90 - (180 / np.pi) * phir

    # Freq. bands for variance contributions
    ig = np.max(np.where(frequencies <= infra_gravity_cutoff))
    # low freq, infragravity, high-freq
    if 1 < ff:
        ublo = np.sqrt(2 * np.sum(Guv[1:ff] * df))
    else:
        ublo = 0
    if ig > ff:
        ubig = np.sqrt(2 * np.sum(Guv[ff:ig] * df))
    else:
        ubig = 0
    if lf < fft_length:
        ubhi = np.sqrt(2 * np.sum(Guv[lf:] * df))
    else:
        ubhi = 0

    ws = {
        'Hrmsp': Hrmsp,
        'Hrmsu': Hrmsu,
        'ubr': ubr,
        'ubr_check': ubr_check,
        'omegar': omegar,
        'Tr': Tr,
        'Tpp': Tpp,
        'Tpu': Tpu,
        'phir': phir,
        'azr': azr,
        'ublo': ublo,
        'ubhi': ubhi,
        'ubig': ubig,
    }

    if check_variances:
        variance_preserved = test_variances(u, v, pressure, Gpp, Guu, Gvv, df, allowable_error=variance_error)
        ws['variance_test_passed'] = variance_preserved

    if show_diagnostic_plot:
        fig, ax = plot_spectra(Guu, Gvv, Guv, Gpp, frequencies, first_frequency_cutoff, ff, last_frequency_cutoff, lf,
                               infra_gravity_cutoff, ig)
        ws['figure'] = fig
        ws['axis'] = ax

    return ws


def plot_spectra(Guu, Gvv, Guv, Gpp, frequencies, first_frequency_cutoff, ff, last_frequency_cutoff, lf,
                 infra_gravity_cutoff, ig):
    """
    internal helper function to plot spectra for diagnostics

        array Guu: east spectra
        array Gvv: north spectra
        array Guv: combined spectra
        array Gpp: pressure spectra
        array frequencies: frequencies
        int ff: first frequency cutoff index
        int lf: last frequency cutoff index
        int ig: infra-gravity wave frequency cutoff index
    :return:  figure and axis handles
    """

    fig, ax = plt.subplots(2, 1, figsize=(15, 5))
    ax[0].plot(frequencies, Guv, label='Guv')
    ax[0].plot(frequencies, Guu, label='Guu')
    ax[0].plot(frequencies, Gvv, label='Gvv')
    ylims = ax[0].get_ylim()
    ax[0].plot([frequencies[ff], frequencies[ff]], ylims, label='ff = %3.1f s' % (1 / frequencies[ff]))
    ax[0].plot([frequencies[lf], frequencies[lf]], ylims, label='lf = %3.1f s' % (1 / frequencies[lf]))
    ax[0].plot([frequencies[ig], frequencies[ig]], ylims, label='ig = %3.1f s' % (1 / frequencies[ig]))
    plt.text(0.5, 0.9,
             'first_frequency_cutoff is {} found at #{} = {} Hz'.format(first_frequency_cutoff, ff, frequencies[ff]),
             transform=ax[0].transAxes)
    plt.text(0.5, 0.8,
             'last_frequency_cutoff is {} found at #{} = {} Hz'.format(last_frequency_cutoff, lf, frequencies[lf]),
             transform=ax[0].transAxes)
    plt.text(0.5, 0.7,
             'infra-gravity cutoff is {} found at #{} = {} Hz'.format(infra_gravity_cutoff, ig, frequencies[ig]),
             transform=ax[0].transAxes)
    ax[0].legend()
    ax[1].plot(frequencies, Gpp, label='Gpp')
    ylims = ax[1].get_ylim()
    ax[1].plot([frequencies[ff], frequencies[ff]], ylims, label='ff = %3.1f s' % (1 / frequencies[ff]))
    ax[1].plot([frequencies[lf], frequencies[lf]], ylims, label='lf = %3.1f s' % (1 / frequencies[lf]))
    ax[1].plot([frequencies[ig], frequencies[ig]], ylims, label='ig = %3.1f s' % (1 / frequencies[ig]))
    ax[1].legend()
    plt.show()

    return fig, ax


def test_variances(u, v, p, Gpp, Guu, Gvv, df, allowable_error=0.0):
    """
    helper function to test for variance preservation
    :param array_like u: velocity time series
    :param array_like v: velocity time series
    :param array_like p: pressure time series
    :param array_like Gpp: pressure spectra
    :param array_like Guu: veloctiy spectra
    :param array_like Gvv: velocity spectra
    :param float df: frequency resolution
    :param float allowable_error:% of allowable difference between the time domain variance and the spectral energy
    :return: bool, True if error limit is met

    """
    result = False

    print('These pairs of time-domain and spectral stats should be very close to equal:')
    varp = np.var(p)
    varP = np.sum(Gpp*df)
    print('var(p) {} == sum(Gpp*df) {}'.format(varp, varP))
    if np.abs((varp-varP)/varp * 100) > allowable_error:
        result = False
    varu = np.var(u)
    varU = np.sum(Guu*df)
    print('var(u) {} == sum(Guu*df) {}'.format(varu, varU))
    if np.abs((varu-varU)/varu * 100) > allowable_error:
        result = False
    varv = np.var(v)
    varV = np.sum(Gvv*df)
    print('var(v) {} == sum(Gvv*df) {}'.format(varv, varV))
    if np.abs((varv-varV)/varv * 100) > allowable_error:
        result = False
    print('sqrt(varv) {} == sqrt(mean(v^2) {}'.format(np.sqrt(varv), np.sqrt(np.mean(v**2))))
    varuv = varu + varv
    varUV = varU + varV
    print('sqrt(2.*(varu+varv)) {} == '.format(np.sqrt(2.*varuv)) +
          'np.sqrt(2*(np.sum(Guu*df) + np.sum(Gvv*df))) {}'.format(np.sqrt(2.*varUV)))
    percent_error = np.abs((np.sqrt(2.*varuv)-np.sqrt(2.*varUV))/np.sqrt(2.*varuv) * 100)
    print(f'percent_error = {percent_error}')
    print(f'allowable_error = {allowable_error}')
    # careful here, it seems that one can't test for 0.0 == 0.0 via the variables
    if (allowable_error == 0.0) & (percent_error == 0.0):
        result = True
    elif (allowable_error > 0.0) & (percent_error < allowable_error):
        result = True

    return result
