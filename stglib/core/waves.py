import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.signal as spsig
import xarray as xr
from tqdm import tqdm

from ..lib import pyDIWASP


def make_diwasp_inputs(
    ds,
    data_type=None,
    layout=None,
    method="IMLM",
    dres=180,
    nsegs=16,
    iter=50,
    freqs=None,
    nfft=None,
    xdir=90,
    dunit="naut",
    nsamps=None,
    smooth="ON",
):
    """Create inputs need to run pyDIWASP processing"""

    # check for user options in ds.attrs
    if "diwasp_method" in ds.attrs:
        method = ds.attrs["diwasp_method"]

    if "diwasp_nsegs" in ds.attrs:
        nsegs = int(ds.attrs["diwasp_nsegs"])

    if "diwasp_dres" in ds.attrs:
        dres = int(ds.attrs["diwasp_dres"])

    if "diwasp_iter" in ds.attrs:
        iter = int(ds.attrs["diwasp_iter"])

    if "diwasp_xdir" in ds.attrs:
        xdir = ds.attrs["diwasp_xdir"]

    if "diwasp_dunit" in ds.attrs:
        dunit = ds.attrs["diwasp_dunit"]

    if "diwasp_nfft" in ds.attrs:
        nfft = int(ds.attrs["diwasp_nfft"])

    if "diwasp_smooth" in ds.attrs:
        smooth = ds.attrs["diwasp_smooth"]

    ID = {}
    ID["fs"] = 1 / float(ds.attrs["sample_interval"])

    if "P_1ac" in ds:
        ID["depth"] = float(
            ds["P_1ac"].mean(dim="sample").values
            + ds.attrs["initial_instrument_height"]
        )

    elif "brangeAST" in ds and ds.attrs["orientation"].lower() == "up":
        ID["depth"] = float(
            ds["brangeAST"].mean(dim="sample").values
            + ds.attrs["initial_instrument_height"]
        )

    else:
        ID["depth"] = float(ds.attrs["WATER_DEPTH"])

    if data_type == "puv":
        ID["datatypes"] = ["pres", "velx", "vely"]

    elif data_type == "suv":
        ID["datatypes"] = ["elev", "velx", "vely"]

    elif data_type == "elev":
        ID["datatypes"] = ["elev"]

    elif data_type == "pres":
        ID["datatypes"] = ["pres"]

    if np.any(layout):
        ID["layout"] = layout
    else:
        raise KeyError(
            "Required DIWASP input parameter <layout> has not been explicitly passed to make_diwasp_inputs def"
        )

    nyfreq = float(ID["fs"]) / 2

    if nsamps is None:
        nsamps = ds.attrs["wave_interval"] * ds.attrs["sample_rate"]

    if nfft is None:
        nfft = next_power_of_2(int(nsamps / nsegs))

    if freqs is None:
        nfreqs = nfft / 2
        flo = np.round(
            1 / (nsamps / ID["fs"] / 32.0), 3
        )  # set minimum frequency to length of wave burst samples use divided by 32
        if nyfreq > 2:
            fhi = 2
        else:
            fhi = nyfreq
        freqs = np.arange(flo, fhi, (fhi - flo) / nfreqs)

    SM = {}
    SM["freqs"] = freqs
    SM["dirs"] = np.arange(0, 360, 360 / dres)
    SM["xaxisdir"] = xdir
    SM["dunit"] = dunit

    EP = {}
    EP["method"] = method
    EP["iter"] = iter
    EP["nfft"] = nfft
    EP["dres"] = int(dres)
    EP["smooth"] = smooth

    return ID, SM, EP


def make_diwasp_puv_suv(ds, layout=None, data_type=None, freqs=None, ibin=0):
    """Calculate Directional Wave Statistic using PyDIWASP"""

    # check number of samples to use in each burst for wave stats (default = all, i.e. samples_per_burst)
    samples_per_burst = len(ds["sample"])  # samples_per_burst
    if "diwasp_nsamps" in ds.attrs:
        nsamps = ds.attrs["diwasp_nsamps"]
    elif "diwasp_pow2" in ds.attrs:  # make diwasp use power of 2 nsamps
        if ds.attrs["diwasp_pow2"].lower() == "true":
            nsamps = floor_power_of_2(samples_per_burst)
        else:
            nsamps = samples_per_burst
    else:
        nsamps = samples_per_burst

    times = []
    s = []
    spec = []
    Hs = []
    Tp = []
    DTp = []
    Dp = []
    # Dm = []
    Tm = []

    dwspid = []

    for burst in tqdm(range(len(ds.time))):
        if "P_1ac" in ds:
            p = ds["P_1ac"].isel(time=burst, sample=range(nsamps))
            if p.isnull().any():
                # Fill NaNs in burst data if possible
                p = var_wave_burst_fill_nans(ds, p)
        else:
            p = None

        if "brangeAST" in ds:
            ast = ds["brangeAST"].isel(time=burst, sample=range(nsamps))
            if ast.isnull().any():
                # Fill NaNs in burst data if possible
                ast = var_wave_burst_fill_nans(ds, ast)

        else:
            ast = None

        u = ds["u_1205"].isel(time=burst, z=ibin, sample=range(nsamps))
        if u.isnull().any():
            # Fill NaNs in burst data if possible
            u = var_wave_burst_fill_nans(ds, u)

        v = ds["v_1206"].isel(time=burst, z=ibin, sample=range(nsamps))
        if v.isnull().any():
            # Fill NaNs in burst data if possible
            v = var_wave_burst_fill_nans(ds, v)

        ID, SM, EP = make_diwasp_inputs(
            ds.isel(time=burst),
            layout=layout,
            data_type=data_type,
            freqs=freqs,
            nsamps=nsamps,
        )
        if data_type == "suv":
            if ast is not None:
                ID["data"] = np.array([ast, u, v]).transpose()
            else:
                raise ValueError(
                    f"Acoustic surface tracking (ast) variable not found cannot continue with {data_type} directional wave analysis"
                )

        elif data_type == "puv":
            if p is not None:
                ID["data"] = np.array([p, u, v]).transpose()
            else:
                raise ValueError(
                    f"Corrected pressure (P_1ac) variable not found cannot continue with {data_type} directional wave analysis"
                )

        elif data_type == "optimized":
            """Optimized method uses suv when ast contains no NaNs and puv when it does"""
            # make list of which method is used
            if ast.isnull().any():
                if p is not None:
                    ID["datatypes"] = ["pres", "velx", "vely"]
                    ID["data"] = np.array([p, u, v]).transpose()
                    dwspid.append("puv")
                else:
                    raise ValueError(
                        f"Corrected pressure (P_1ac) variable not found cannot continue with {data_type} directional wave analysis"
                    )

            else:
                if ast is not None:
                    ID["datatypes"] = ["elev", "velx", "vely"]
                    ID["data"] = np.array([ast, u, v]).transpose()
                    dwspid.append("suv")
                else:
                    raise ValueError(
                        f"Acoustic surface tracking (ast) variable not found cannot continue with {data_type} directional wave analysis"
                    )

        opts = ["MESSAGE", 0, "PLOTTYPE", 0]

        # check to make sure no NaNs in input data
        if np.any(np.isnan(ID["data"])):

            SMout = {}
            SMout["S"] = np.full((len(SM["freqs"]), len(SM["dirs"])), np.nan)
            SMout["freqs"] = SM["freqs"]
            SMout["dirs"] = SM["dirs"]
            WVout = np.full(4, np.nan)

        else:

            [SMout, EPout] = pyDIWASP.dirspec.dirspec(ID, SM, EP, opts)
            WVout = pyDIWASP.infospec.infospec(SMout)

        #  Extract real part of directional spectra and calculate frequency spectra
        Dnn = np.real(SMout["S"])
        Snn = np.trapezoid(Dnn, axis=1, x=SMout["dirs"])

        # Make tail if using puv
        if dwspid:
            print(dwspid[-1])
        if data_type == "puv" or (dwspid and dwspid[-1] == "puv"):
            print("Apply cutoff and add tail for puv method")
            if "pressure_sensor_height" in ds.attrs:
                z = ds.attrs["pressure_sensor_height"]
            else:
                warnings.warn(
                    "pressure_sensor_height not specified; using initial_instrument_height to compute wave statistics"
                )
                z = ds.attrs["initial_instrument_height"]

            h = np.mean(p.values) + z
            f = SMout["freqs"]
            k = qkfs(2 * np.pi * f, h)
            Kp = transfer_function(k, h, z)

            # Find cutoff
            tailind, noisecutind, fpeakcutind, Kpcutind = define_cutoff(
                f, Snn * (Kp**2), Kp, noise=0.9
            )

            if ~np.isnan(tailind):
                Snn = make_tail(f, Snn, tailind)

                # Make tail for dspec
                Dnn = make_dspec_tail(f, SMout["dirs"], Dnn, tailind)

        m0 = make_moment(SMout["freqs"], Snn, 0)
        m2 = make_moment(SMout["freqs"], Snn, 2)
        mwp = make_Tm(m0, m2)
        hs = make_Hs(m0)
        tp = make_Tp(
            xr.DataArray(Snn, dims=["frequency"], coords={"frequency": SMout["freqs"]})
        )

        # append outputs to list
        times.append(ds["time"][burst].values)
        s.append(Dnn)
        spec.append(Snn)
        Hs.append(hs)
        Tp.append(tp)
        DTp.append(WVout[2])
        Dp.append(WVout[3])
        Tm.append(mwp)

    dspec = np.stack(s, axis=0)
    fspec = np.stack(spec, axis=0)

    dwv = xr.Dataset()
    dwv["time"] = xr.DataArray(np.array(times), dims="time")
    dwv["time"] = pd.DatetimeIndex(dwv["time"])
    dwv["diwasp_frequency"] = xr.DataArray(SMout["freqs"], dims="diwasp_frequency")
    dwv["diwasp_direction"] = xr.DataArray(SMout["dirs"], dims="diwasp_direction")
    dwv["diwasp_dspec"] = xr.DataArray(
        dspec, dims=["time", "diwasp_frequency", "diwasp_direction"]
    )

    # dwv["diwasp_fspec"] = xr.DataArray(
    #    np.trapezoid(dspec, axis=2, x=SMout["dirs"]),
    #    dims=["time", "diwasp_frequency"],
    # )

    dwv["diwasp_fspec"] = xr.DataArray(
        fspec,
        dims=["time", "diwasp_frequency"],
    )

    # find mean wave direction
    Dm = make_mwd(
        dwv["diwasp_frequency"],
        dwv["diwasp_direction"],
        dwv["diwasp_dspec"],
        diwasp=True,
    )

    dwv["diwasp_hs"] = xr.DataArray(np.array(Hs), dims="time")
    dwv["diwasp_tp"] = xr.DataArray(np.array(Tp), dims="time")
    dwv["diwasp_dtp"] = xr.DataArray(np.array(DTp), dims="time")
    dwv["diwasp_tm"] = xr.DataArray(np.array(Tm), dims="time")
    dwv["diwasp_dp"] = xr.DataArray(np.array(Dp), dims="time")
    dwv["diwasp_dm"] = xr.DataArray(np.round(Dm.values, 0), dims="time")

    # Add List of input data type for optimized method
    print(np.array(dwspid))
    if dwspid:
        dwv["diwasp_type"] = xr.DataArray(np.array(dwspid), dims="time")

    # make some diwasp attrs
    if "diwasp_ibin" not in ds.attrs:
        dwv.attrs["diwasp_ibin"] = ibin

    # if "diwasp_inputs" not in ds.attrs:
    if data_type == "optimized":
        dwv.attrs["diwasp_inputs"] = (
            "optimized for ['elev', 'velx', 'vely'] or ['pres'. 'velx', 'vely']"
        )
    else:
        dwv.attrs["diwasp_inputs"] = ID["datatypes"]

    if "diwasp_method" not in ds.attrs:
        dwv.attrs["diwasp_method"] = EP["method"]
    if "diwasp_nsamps" not in ds.attrs:
        dwv.attrs["diwasp_nsamps"] = nsamps
    if "diwasp_nfft" not in ds.attrs:
        dwv.attrs["diwasp_nfft"] = EP["nfft"]
    if "diwasp_dres" not in ds.attrs:
        dwv.attrs["diwasp_dres"] = EP["dres"]
    if "diwasp_iter" not in ds.attrs:
        dwv.attrs["diwasp_iter"] = EP["iter"]
    if "diwasp_xdir" not in ds.attrs:
        dwv.attrs["diwasp_xdir"] = SM["xaxisdir"]
    if "diwasp_dunit" not in ds.attrs:
        dwv.attrs["diwasp_dunit"] = SM["dunit"]

    return dwv


def make_diwasp_elev_pres(ds, layout=None, data_type=None, freqs=None):
    """Calculate Directional Wave Statistic using PyDIWASP"""

    # check number of samples to use in each burst for wave stats (default = all, i.e. samples_per_burst)
    samples_per_burst = len(ds["sample"])  # samples_per_burst
    if "diwasp_nsamps" in ds.attrs:
        nsamps = ds.attrs["diwasp_nsamps"]
    elif "diwasp_pow2" in ds.attrs:  # make diwasp use power of 2 nsamps
        if ds.attrs["diwasp_pow2"].lower() == "true":
            nsamps = floor_power_of_2(samples_per_burst)
        else:
            nsamps = samples_per_burst
    else:
        nsamps = samples_per_burst

    times = []
    spec = []
    Hs = []
    Tp = []
    Tm = []

    dwspid = []

    for burst in tqdm(range(len(ds.time))):

        if "P_1ac" in ds:
            p = ds["P_1ac"].isel(time=burst, sample=range(nsamps))
            if p.isnull().any():
                # Fill NaNs in burst data if possible
                p = var_wave_burst_fill_nans(ds, p)
        else:
            p = None

        if "brangeAST" in ds:
            ast = ds["brangeAST"].isel(time=burst, sample=range(nsamps))
            if ast.isnull().any():
                # Fill NaNs in burst data if possible
                ast = var_wave_burst_fill_nans(ds, ast)
        else:
            ast = None

        ID, SM, EP = make_diwasp_inputs(
            ds.isel(time=burst),
            data_type=data_type,
            layout=layout,
            freqs=freqs,
            nsamps=nsamps,
        )
        if data_type == "elev":
            if ast is not None:
                ID["data"] = np.atleast_2d(ast).transpose()
            else:
                raise ValueError(
                    f"Acoustic surface tracking (ast) variable not found cannot continue with {data_type} wave analysis"
                )

        elif data_type == "pres":
            if p is not None:
                ID["data"] = np.atleast_2d(p).transpose()
            else:
                raise ValueError(
                    f"Correct pressure (P_1ac) variable not found cannot continue with {data_type} wave analysis"
                )

        elif data_type == "optimized-nd":
            """Optimized-ND method uses elev when ast contains no NaNs and pres when it does"""
            # make list of which method is used
            if ast.isnull().any():
                if p is not None:
                    ID["datatypes"] = ["pres"]
                    ID["data"] = np.atleast_2d(p).transpose()
                    dwspid.append("pres")
                else:
                    raise ValueError(
                        f"Corrected pressure (P_1ac) variable not found cannot continue with {data_type} wave analysis"
                    )

            else:
                if ast is not None:
                    ID["datatypes"] = ["elev"]
                    ID["data"] = np.atleast_2d(ast).transpose()
                    dwspid.append("elev")
                else:
                    raise ValueError(
                        f"Acoustic surface tracking (ast) variable not found cannot continue with {data_type} wave analysis"
                    )

        opts = ["MESSAGE", 0, "PLOTTYPE", 0]

        # check to make sure no NaNs in input data
        if np.any(np.isnan(ID["data"])):

            SMout = {}
            SMout["S"] = np.full((len(SM["freqs"]), len(SM["dirs"])), np.nan)
            SMout["freqs"] = SM["freqs"]
            SMout["dirs"] = SM["dirs"]
            # WVout = np.full(4, np.nan)

        else:

            [SMout, EPout] = pyDIWASP.dirspec.dirspec(ID, SM, EP, opts)
            # WVout = pyDIWASP.infospec.infospec(SMout)

        #  Extract real part of directional spectra and calculate frequency spectra
        Dnn = np.real(SMout["S"])
        Snn = np.trapezoid(Dnn, axis=1, x=SMout["dirs"])

        # Make tail if using pres
        if data_type == "pres" or (dwspid and dwspid[-1] == "pres"):
            print("Apply cutofff and add tail for pres method")
            if "pressure_sensor_height" in ds.attrs:
                z = ds.attrs["pressure_sensor_height"]
            else:
                warnings.warn(
                    "pressure_sensor_height not specified; using initial_instrument_height to compute wave statistics"
                )
                z = ds.attrs["initial_instrument_height"]

            h = np.mean(p.values) + z
            f = SMout["freqs"]
            k = qkfs(2 * np.pi * f, h)
            Kp = transfer_function(k, h, z)

            # Find cutoff
            tailind, noisecutind, fpeakcutind, Kpcutind = define_cutoff(
                f, Snn * (Kp**2), Kp, noise=0.9
            )

            if ~np.isnan(tailind):
                Snn = make_tail(f, Snn, tailind)

        # calculate mean period
        m0 = make_moment(SMout["freqs"], Snn, 0)
        m2 = make_moment(SMout["freqs"], Snn, 2)
        mwp = make_Tm(m0, m2)

        hs = make_Hs(m0)
        tp = make_Tp(
            xr.DataArray(Snn, dims=["frequency"], coords={"frequency": SMout["freqs"]})
        )

        # append outputs to list
        times.append(ds["time"][burst].values)
        spec.append(Snn)
        Hs.append(hs)
        Tp.append(tp)
        Tm.append(mwp)

    fspec = np.stack(spec, axis=0)

    dwv = xr.Dataset()
    dwv["time"] = xr.DataArray(np.array(times), dims="time")
    dwv["time"] = pd.DatetimeIndex(dwv["time"])
    dwv["diwasp_frequency"] = xr.DataArray(SMout["freqs"], dims="diwasp_frequency")
    dwv["diwasp_fspec"] = xr.DataArray(
        fspec,
        dims=["time", "diwasp_frequency"],
    )

    dwv["diwasp_hs"] = xr.DataArray(np.array(Hs), dims="time")
    dwv["diwasp_tp"] = xr.DataArray(np.array(Tp), dims="time")
    dwv["diwasp_tm"] = xr.DataArray(np.array(Tm), dims="time")

    # Add List of input data type for optimized-nd method
    print(np.array(dwspid))
    if dwspid:
        dwv["diwasp_type"] = xr.DataArray(np.array(dwspid), dims="time")

    # if "diwasp_inputs" not in ds.attrs:
    if data_type == "optimized-nd":
        dwv.attrs["diwasp_inputs"] = "optimized for ['elev'] or ['pres']"
    else:
        dwv.attrs["diwasp_inputs"] = ID["datatypes"]

    if "diwasp_method" not in ds.attrs:
        dwv.attrs["diwasp_method"] = EP["method"]
    if "diwasp_nsamps" not in ds.attrs:
        dwv.attrs["diwasp_nsamps"] = nsamps
    if "diwasp_nfft" not in ds.attrs:
        dwv.attrs["diwasp_nfft"] = EP["nfft"]
    if "diwasp_dres" not in ds.attrs:
        dwv.attrs["diwasp_dres"] = EP["dres"]
    if "diwasp_iter" not in ds.attrs:
        dwv.attrs["diwasp_iter"] = EP["iter"]
    if "diwasp_xdir" not in ds.attrs:
        dwv.attrs["diwasp_xdir"] = SM["xaxisdir"]
    if "diwasp_dunit" not in ds.attrs:
        dwv.attrs["diwasp_dunit"] = SM["dunit"]

    return dwv


def make_waves_ds(ds):
    print("Computing wave statistics")
    if "P_1ac" in ds:
        presvar = "P_1ac"
    else:
        warnings.warn(
            "atmospherically corrected pressure not available; using raw pressure to compute wave statistics"
        )
        presvar = "P_1"

    if "spec_nsegs" in ds.attrs:
        nsegs = ds.attrs["spec_nsegs"]
    else:
        # set to default = 16
        nsegs = 16

    nsamps = len(ds["sample"].values)
    nfft = next_power_of_2(int(nsamps / nsegs))

    f, Pxx = pressure_spectra(
        ds[presvar].squeeze().values,
        fs=1 / ds.attrs["sample_interval"],
        nperseg=nfft,
    )

    if "pressure_sensor_height" in ds.attrs:
        z = ds.attrs["pressure_sensor_height"]
    else:
        warnings.warn(
            "pressure_sensor_height not specified; using initial_instrument_height to compute wave statistics"
        )
        z = ds.attrs["initial_instrument_height"]
    h = ds[presvar].squeeze().mean(dim="sample") + z

    k = np.asarray([qkfs(2 * np.pi * f, x) for x in h.values])

    Kp = transfer_function(k, h, z)
    Pnn = elevation_spectra(Pxx, Kp)

    spec = xr.Dataset()

    spec["Pnn"] = xr.DataArray(Pnn, dims=("time", "frequency"), coords=(ds["time"], f))
    spec["Pxx"] = xr.DataArray(Pxx, dims=("time", "frequency"), coords=(ds["time"], f))
    spec["Kp"] = xr.DataArray(Kp, dims=("time", "frequency"), coords=(ds["time"], f))
    spec["k"] = xr.DataArray(k, dims=("time", "frequency"), coords=(ds["time"], f))

    if "waves_fractional_noise" in ds.attrs:
        noise = ds.attrs["waves_fractional_noise"]
    else:
        noise = 0.9
    tailind, noisecutind, fpeakcutind, Kpcutind = zip(
        *[
            define_cutoff(f, x, y, noise=noise)
            for x, y in zip(spec["Pxx"].values, spec["Kp"].values)
        ]
    )
    spec["tailind"] = xr.DataArray(np.asarray(tailind), dims="time")
    spec["noisecutind"] = xr.DataArray(np.asarray(noisecutind), dims="time")
    spec["fpeakcutind"] = xr.DataArray(np.asarray(fpeakcutind), dims="time")
    spec["Kpcutind"] = xr.DataArray(np.asarray(Kpcutind), dims="time")
    thetail = [
        make_tail(
            spec["frequency"],
            spec["Pnn"].sel(time=time),
            spec["tailind"].sel(time=time).values,
        )
        for time in spec["time"]
    ]
    spec["pspec"] = xr.DataArray(thetail, dims=("time", "frequency"))

    # clip wave spectral frequencies if user enabled - must be done before taking any spectral moments
    if "clip_f" in ds.attrs and ds.attrs["clip_f"].lower() == "true":
        spec = clip_f(spec, nsamps=nsamps, fs=1 / ds.attrs["sample_interval"])

    spec["m0"] = xr.DataArray(
        make_moment(spec["frequency"], spec["pspec"], 0), dims="time"
    )
    spec["m2"] = xr.DataArray(
        make_moment(spec["frequency"], spec["pspec"], 2), dims="time"
    )
    spec["wh_4061"] = xr.DataArray(make_Hs(spec["m0"]), dims="time")
    spec["wp_4060"] = xr.DataArray(make_Tm(spec["m0"], spec["m2"]), dims="time")
    spec["wp_peak"] = xr.DataArray(make_Tp(spec["pspec"]), dims="time")

    return spec


def make_waves_ds_elev(ds):
    """Calculate wave statistics using sea-surface elevation"""
    print("Computing wave statistics")
    if "elev" in ds:
        var = "elev"
    elif "brange" in ds:
        var = "brange"
    elif "brangeAST" in ds:
        var = "brangeAST"

    if "spec_nsegs" in ds.attrs:
        nsegs = ds.attrs["spec_nsegs"]
    else:
        # set to default = 16
        nsegs = 16

    nsamps = len(ds["sample"].values)
    nfft = next_power_of_2(int(nsamps / nsegs))

    f, Pxx = pressure_spectra(
        ds[var].squeeze().values,
        fs=1 / ds.attrs["sample_interval"],
        nperseg=nfft,
    )

    # trim frequencies
    nyfreq = float(ds.attrs["sample_rate"]) / 2

    # set minimum frequency so that 32 complete intervals are contained in the wave burst
    # set maximum frequency to lesser of Nyquist (sample_rate/2) or 2Hz
    # These frequency settings maximize efficiency of frequency space over the resolvable range of the measurements.
    flo = np.round(1 / (nsamps / ds.attrs["sample_rate"] / 32.0), 3)
    if nyfreq > 2:
        fhi = 2
    else:
        fhi = nyfreq

    ind = (f >= flo) & (f <= fhi)
    f = f[ind]
    Pxx = Pxx[:, ind]

    Pnn = Pxx  # measurement is of the sea-surface directly

    spec = xr.Dataset()

    spec["Pnn"] = xr.DataArray(Pnn, dims=("time", "frequency"), coords=(ds["time"], f))

    # Use Pnn directly, no cutoff and no tail
    spec["sspec"] = xr.DataArray(Pnn, dims=("time", "frequency"))
    spec["m0"] = xr.DataArray(
        make_moment(spec["frequency"], spec["sspec"], 0), dims="time"
    )
    spec["m2"] = xr.DataArray(
        make_moment(spec["frequency"], spec["sspec"], 2), dims="time"
    )
    spec["wh_4061"] = xr.DataArray(make_Hs(spec["m0"]), dims="time")
    spec["wp_4060"] = xr.DataArray(make_Tm(spec["m0"], spec["m2"]), dims="time")
    spec["wp_peak"] = xr.DataArray(make_Tp(spec["sspec"]), dims="time")

    return spec


def pressure_spectra(x, fs=1.0, window="hann", nperseg=256, **kwargs):
    """Compute pressure spectral density using Welch's method

    Parameters
    ----------
    x : array_like
        Time-series of pressure data
    fs : float, optional
        Sampling frequency (Hz)
    window : str, optional
        Window, default 'hann'
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
    """Compute elevation spectra using linear wave theory and transfer function"""
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
        Pressure transfer function
    """
    if isinstance(h, float) or h.ndim == 0:
        Kp = np.cosh(k * z) / np.cosh(k * h)
    else:
        Kp = np.cosh(k * z) / np.cosh(k * np.expand_dims(h, 1))

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

    noisecut = 12 * np.mean(Pxx[f >= noise * f[-1]])
    # Look for above cut-off freqs and take highest
    tmp = np.where(Pxx > noisecut)[0]
    if len(tmp) == 0:
        noisecutind = 0
    else:
        noisecutind = tmp[-1]

    fpeakcut = 1.1 * f[np.argmax(Pxx)]
    fpeakcutind = np.searchsorted(f, fpeakcut)  # cutoff based on 1.1*fp
    try:
        Kpcutind = np.nonzero(Kp > 0.1)[0][-1] + 1  # cutoff keeping only Kp > 0.1
    except IndexError:
        Kpcutind = 0

    # take the more conservative of either K_p or noise cutoff
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
        tail[ti:] = Pnn[ti] * (f[ti:] / f[ti]) ** -4
        return np.hstack((Pnn[:ti], tail[ti:]))


def make_dspec_tail(f, dirs, dspec, tailind):
    """Add tail to directional spectra"""

    dspectail = []
    for k in range(len(dirs)):
        Pnn = dspec[:, k]
        thetail = make_tail(f, Pnn, tailind)
        dspectail.append(thetail)

    return np.stack(dspectail, axis=1)


def make_mwd(freqs, dirs, dspec, diwasp=False):
    """Create mean wave direction (EPIC 4062) variable"""

    Sxsin = dspec * np.sin(np.deg2rad(dirs))
    Sxcos = dspec * np.cos(np.deg2rad(dirs))

    if diwasp:
        Dnum = Sxsin.integrate("diwasp_frequency").integrate("diwasp_direction")
        Ddnom = Sxcos.integrate("diwasp_frequency").integrate("diwasp_direction")
    else:
        Dnum = Sxsin.integrate("frequency").integrate("direction")
        Ddnom = Sxcos.integrate("frequency").integrate("direction")

    Dm = np.rad2deg(np.arctan(np.abs(Dnum / Ddnom)))

    Dm[(Dnum > 0) & (Ddnom < 0)] = 180 - Dm[(Dnum > 0) & (Ddnom < 0)]

    Dm[(Dnum < 0) & (Ddnom < 0)] = 180 + Dm[(Dnum < 0) & (Ddnom < 0)]

    Dm[(Dnum < 0) & (Ddnom > 0)] = 360 - Dm[(Dnum < 0) & (Ddnom > 0)]

    return Dm


def make_moment(f, Pnn, n):
    """Compute nth moment (m0, m1, m2, etc.) of power spectra"""
    return np.trapezoid(Pnn * f**n, x=f)


def make_Hs(m0):
    """Compute significant wave height as 4 * np.sqrt(m0)"""
    return 4 * np.sqrt(m0)


def make_Tm(m0, m2):
    """Compute mean period as sqrt(m0 / m2)"""
    return np.sqrt(m0 / m2)


def make_Tp(Pnn):
    """Compute peak period as 1 / fp,
    where fp is the frequency with greatest energy in the elevation spectra"""
    # ensure we don't return 0 frequency as a peak period
    fp = Pnn["frequency"][Pnn.fillna(0).argmax(dim="frequency")].values
    fp[fp == 0] = np.nan
    return 1 / fp


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
    y = y - ((y * t - x) / (t + y * (1 - t**2)))
    t = np.tanh(y)
    y = y - ((y * t - x) / (t + y * (1 - t**2)))
    t = np.tanh(y)
    y = y - ((y * t - x) / (t + y * (1 - t**2)))

    return y / h


def detrend_nan(arr):
    """Helper function for spsig.detrend().
    This fills nans with zeros so that detrend won't fail, and then nans out the
    values that had nans in them (since the result of the detrend is not reliable
    and we won't be using the values anyway)"""

    goods = ~arr.isnull()
    arr = arr.where(goods, 0)
    arr = spsig.detrend(arr)
    # now it's a numpy array so can't use .where
    arr[~goods] = np.nan
    return arr


def puv_quick_vectorized(
    pressure,
    u,
    v,
    depth,
    height_of_pressure,
    height_of_velocity,
    sampling_frequency,
    fft_length=512,
    rho=1025.0,
    first_frequency_cutoff=1 / 50,
    infra_gravity_cutoff=0.05,
    last_frequency_cutoff=1 / 5,
    fft_window_type="hann",
    show_diagnostic_plot=False,
    check_variances=False,
    variance_error=0.0,
    overlap_length="default",
):
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

    Daniel Nowacki vectorized version for speed improvements
    """

    gravity = 9.81  # m/s^2
    if fft_window_type == "hann":
        fft_window_type = "hann"  # this is just the way scipy signal likes it
    if overlap_length == "default":
        overlap_length = int(np.floor(fft_length / 2))

    pressure = detrend_nan(pressure)
    u = detrend_nan(u)
    v = detrend_nan(v)

    # compute wave height from velocities

    # Determine velocity spectra for u and v
    [frequencies, Gpp] = spsig.welch(
        rho * gravity * pressure,
        fs=sampling_frequency,
        window=fft_window_type,
        nperseg=fft_length,
        noverlap=overlap_length,
    )

    df = frequencies[2] - frequencies[1]
    [_, Guu] = spsig.welch(
        u,
        fs=sampling_frequency,
        window=fft_window_type,
        nperseg=fft_length,
        noverlap=overlap_length,
    )
    [_, Gvv] = spsig.welch(
        v,
        fs=sampling_frequency,
        window=fft_window_type,
        nperseg=fft_length,
        noverlap=overlap_length,
    )

    # determine wave number
    omega = np.array(
        [2 * np.pi * x for x in frequencies]
    )  # omega must be numpy array for qkfs
    # catch numpy errors
    np.seterr(divide="ignore", invalid="ignore")
    # k = qkfs(omega, float(depth))  # make sure it is float, or qkfs will bomb
    k = np.array([qkfs(omega, d) for d in depth])
    np.seterr(divide=None, invalid=None)

    # compute linear wave transfer function
    kh = k * np.broadcast_to(np.atleast_2d(depth).T, k.shape)
    kzp = k * height_of_pressure
    kzuv = k * height_of_velocity
    Hp = np.ones(kh.shape)
    Huv = np.ones(kh.shape)

    # change wavenumber at 0 Hz to 1 to avoid divide by zero
    # for some reason in the MATLAB version CRS tests omega for nans instead of k.
    # Here we test k also because that's where the nans show up

    Hp = rho * gravity * (np.cosh(kzp) / np.cosh(kh))
    Huv = omega * (np.cosh(kzuv) / np.sinh(kh))

    # do this after so we can vectorize the math above
    if np.isnan(omega[0]) or (omega[0] <= 0):
        Hp[:, 0] = 1
        Huv[:, 0] = 1

    Hp[np.isnan(k[:, 0]), 0] = 1
    Huv[np.isnan(k[:, 0]), 0] = 1

    # combine horizontal velocity spectra
    Guv = Guu + Gvv

    # create cut off frequency, so noise is not magnified
    # at least in first testing, subtracting 1 here got closer to the intended freq. cutoff value
    ff = np.argmax(frequencies > first_frequency_cutoff) - 1
    lf = np.argmax(frequencies > last_frequency_cutoff)

    # Determine wave height for velocity spectra
    Snp = Gpp[:, ff:lf] / (Hp[:, ff:lf] ** 2)
    Snu = Guv[:, ff:lf] / (Huv[:, ff:lf] ** 2)
    fclip = frequencies[ff:lf]

    Kp = transfer_function(k, depth, height_of_pressure)
    tailind, noisecutind, fpeakcutind, Kpcutind = zip(
        *[define_cutoff(frequencies, x, y) for x, y in zip(Gpp, Kp)]
    )
    tailind = np.array(tailind)
    Snp_tail = [
        make_tail(frequencies, Gpp[burst, :] / Hp[burst, :] ** 2, tailind[burst])
        for burst in range(Gpp.shape[0])
    ]
    Snp_tail = np.array(Snp_tail)
    # skip the zero frequency -- energy persists...
    Snp_tail[:, 0] = np.nan

    Kp_u = transfer_function(k, depth, height_of_velocity)
    tailind_u, noisecutind_u, fpeakcutind_u, Kpcutind_u = zip(
        *[define_cutoff(frequencies, x, y) for x, y in zip(Guv, Kp_u)]
    )
    tailind_u = np.array(tailind_u)
    Snu_tail = [
        make_tail(frequencies, Guv[burst, :] / Huv[burst, :] ** 2, tailind_u[burst])
        for burst in range(Guv.shape[0])
    ]
    Snu_tail = np.array(Snu_tail)
    # skip the zero frequency -- energy persists...
    Snu_tail[:, 0] = np.nan

    # Determine rms wave height (multiply by another sqrt(2) for Hs)
    # Thornton and Guza say Hrms = sqrt(8 mo)
    Hrmsu = 2 * np.sqrt(2 * np.sum(Snu * df, axis=-1))
    Hrmsp = 2 * np.sqrt(2 * np.sum(Snp * df, axis=-1))

    # # skip the zero frequency by starting at 1
    Hrmsu_tail = 2 * np.sqrt(2 * np.sum(Snu_tail[:, 1:] * df, axis=-1))
    Hrmsp_tail = 2 * np.sqrt(2 * np.sum(Snp_tail[:, 1:] * df, axis=-1))

    # These are representative orbital velocities for w-c calculations,
    # according to Madsen (1994) Coastal Engineering 1994, Proc., 24th
    # Intl. Conf., Coastal Eng. Res. Council / ASCE. pp.384-398.
    # (esp. p. 395)
    ubr = np.sqrt(2 * np.sum(Guv[:, ff:lf] * df, axis=-1))
    ubr_check = np.sqrt(2 * np.var(u) + 2 * np.var(v))
    omegar = np.sum(omega[ff:lf] * Guv[:, ff:lf] * df, axis=-1) / np.sum(
        Guv[:, ff:lf] * df, axis=-1
    )
    Tr = 2 * np.pi / omegar

    if len(np.where(np.isnan(Snp))) > 0 | len(np.where(Snp == 0)) > 0:
        Tpp = np.nan
    else:
        jpeak = np.argmax(Snp, axis=-1)  # index location of the maximum value
        Tpp = 1 / fclip[jpeak]

    if len(np.where(np.isnan(Snu))) > 0 | len(np.where(Snu == 0)) > 0:
        Tpu = np.nan
    else:
        jpeak = np.argmax(Snu, axis=-1)
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
    rr = np.corrcoef(u, v).diagonal(u.shape[0])
    ortest = np.sign(rr)
    phir = np.arctan2(
        ortest * np.sum(Gvv[:, ff:lf] * df, axis=-1),
        np.sum(Guu[:, ff:lf] * df, axis=-1),
    )
    phir_tail = np.arctan2(
        ortest * np.sum(Gvv * df, axis=-1), np.sum(Guu * df, axis=-1)
    )

    # convert to degrees; convert to geographic azimuth (0-360, 0=north)
    azr = 90 - (180 / np.pi) * phir
    azr_tail = 90 - (180 / np.pi) * phir_tail

    # Freq. bands for variance contributions
    ig = np.max(np.where(frequencies <= infra_gravity_cutoff))
    # low freq, infragravity, high-freq
    if 1 < ff:
        ublo = np.sqrt(2 * np.sum(Guv[:, 1:ff] * df, axis=-1))
    else:
        ublo = np.zeros_like(Hrmsp)
    if ig > ff:
        ubig = np.sqrt(2 * np.sum(Guv[:, ff:ig] * df, axis=-1))
    else:
        ubig = np.zeros_like(Hrmsp)
    if lf < fft_length:
        ubhi = np.sqrt(2 * np.sum(Guv[:, lf:] * df, axis=-1))
    else:
        ubhi = np.zeros_like(Hrmsp)

    ws = {
        "Hrmsp": Hrmsp,
        "Hrmsu": Hrmsu,
        "ubr": ubr,
        "ubr_check": ubr_check,
        "omegar": omegar,
        "Tr": Tr,
        "Tpp": Tpp,
        "Tpu": Tpu,
        "phir": phir,
        "azr": azr,
        "ublo": ublo,
        "ubhi": ubhi,
        "ubig": ubig,
        "frequencies": frequencies,
        "Gpp": Gpp,
        "Guv": Guv,
        "Guu": Guu,
        "Gvv": Gvv,
        "Snp": Snp,
        "Snu": Snu,
        "Snp_tail": Snp_tail,
        "Snu_tail": Snu_tail,
        "Hrmsp_tail": Hrmsp_tail,
        "Hrmsu_tail": Hrmsu_tail,
        "phir_tail": phir_tail,
        "azr_tail": azr_tail,
        "fclip": fclip,
    }

    if check_variances:
        variance_preserved = test_variances(
            u, v, pressure, Gpp, Guu, Gvv, df, allowable_error=variance_error
        )
        ws["variance_test_passed"] = variance_preserved

    if show_diagnostic_plot:
        fig, ax = plot_spectra(
            Guu,
            Gvv,
            Guv,
            Gpp,
            frequencies,
            first_frequency_cutoff,
            ff,
            last_frequency_cutoff,
            lf,
            infra_gravity_cutoff,
            ig,
        )
        ws["figure"] = fig
        ws["axis"] = ax

    return ws


def puv_quick(
    pressure,
    u,
    v,
    depth,
    height_of_pressure,
    height_of_velocity,
    sampling_frequency,
    fft_length=512,
    rho=1025.0,
    first_frequency_cutoff=1 / 50,
    infra_gravity_cutoff=0.05,
    last_frequency_cutoff=1 / 5,
    fft_window_type="hann",
    show_diagnostic_plot=False,
    check_variances=False,
    variance_error=0.0,
    overlap_length="default",
):
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
    if fft_window_type == "hann":
        fft_window_type = "hann"  # this is just the way scipy signal likes it
    if overlap_length == "default":
        overlap_length = int(np.floor(fft_length / 2))

    pressure = spsig.detrend(pressure)
    u = spsig.detrend(u)
    v = spsig.detrend(v)

    # compute wave height from velocities

    # Determine velocity spectra for u and v
    [frequencies, Gpp] = spsig.welch(
        rho * gravity * pressure,
        fs=sampling_frequency,
        window=fft_window_type,
        nperseg=fft_length,
        noverlap=overlap_length,
    )

    df = frequencies[2] - frequencies[1]
    [_, Guu] = spsig.welch(
        u,
        fs=sampling_frequency,
        window=fft_window_type,
        nperseg=fft_length,
        noverlap=overlap_length,
    )
    [_, Gvv] = spsig.welch(
        v,
        fs=sampling_frequency,
        window=fft_window_type,
        nperseg=fft_length,
        noverlap=overlap_length,
    )

    # determine wave number
    omega = np.array(
        [2 * np.pi * x for x in frequencies]
    )  # omega must be numpy array for qkfs
    # catch numpy errors
    np.seterr(divide="ignore", invalid="ignore")
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
    i = np.array(
        range(nf)
    )  # this is an index, thus needs to start at first element, in this case 0
    # for some reason in the MATLAB version CRS tests omega for nans instead of k.
    # Here we test k also because that's where the nans show up
    if (
        np.isnan(omega[0]) or np.isnan(k[0]) or (omega[0] <= 0)
    ):  # 0 Hz is the first element
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

    Kp = transfer_function(k, depth, height_of_pressure)
    tailind, noisecutind, fpeakcutind, Kpcutind = define_cutoff(frequencies, Gpp, Kp)
    if np.isnan(tailind):
        Snp_tail = np.full_like(frequencies, np.nan)
    else:
        Snp_tail = make_tail(frequencies, Gpp / Hp**2, tailind)
    # skip the zero frequency -- energy persists...
    Snp_tail[0] = np.nan

    Kp_u = transfer_function(k, depth, height_of_velocity)
    tailind_u, noisecutind_u, fpeakcutind_u, Kpcutind_u = define_cutoff(
        frequencies, Guv, Kp_u
    )

    if np.isnan(tailind_u):
        Snu_tail = np.full_like(frequencies, np.nan)
    else:
        Snu_tail = make_tail(frequencies, Guv / Huv**2, tailind_u)
    # skip the zero frequency -- energy persists...
    Snu_tail[0] = np.nan

    # Determine rms wave height (multiply by another sqrt(2) for Hs)
    # Thornton and Guza say Hrms = sqrt(8 mo)
    Hrmsu = 2 * np.sqrt(2 * np.sum(Snu * df))
    Hrmsp = 2 * np.sqrt(2 * np.sum(Snp * df))

    # skip the zero frequency by starting at 1
    if np.isnan(tailind_u):
        Hrmsu_tail = np.nan
    else:
        Hrmsu_tail = 2 * np.sqrt(2 * np.sum(Snu_tail[1:] * df))
    if np.isnan(tailind):
        Hrmsp_tail = np.nan
    else:
        Hrmsp_tail = 2 * np.sqrt(2 * np.sum(Snp_tail[1:] * df))

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
    phir_tail = np.arctan2(ortest * np.sum(Gvv * df), np.sum(Guu * df))

    # convert to degrees; convert to geographic azimuth (0-360, 0=north)
    azr = 90 - (180 / np.pi) * phir
    azr_tail = 90 - (180 / np.pi) * phir_tail

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
        "Hrmsp": Hrmsp,
        "Hrmsu": Hrmsu,
        "ubr": ubr,
        "ubr_check": ubr_check,
        "omegar": omegar,
        "Tr": Tr,
        "Tpp": Tpp,
        "Tpu": Tpu,
        "phir": phir,
        "azr": azr,
        "ublo": ublo,
        "ubhi": ubhi,
        "ubig": ubig,
        "frequencies": frequencies,
        "Gpp": Gpp,
        "Guv": Guv,
        "Guu": Guu,
        "Gvv": Gvv,
        "Snp": Snp,
        "Snu": Snu,
        "Snp_tail": Snp_tail,
        "Snu_tail": Snu_tail,
        "Hrmsp_tail": Hrmsp_tail,
        "Hrmsu_tail": Hrmsu_tail,
        "phir_tail": phir_tail,
        "azr_tail": azr_tail,
        "fclip": fclip,
    }

    if check_variances:
        variance_preserved = test_variances(
            u, v, pressure, Gpp, Guu, Gvv, df, allowable_error=variance_error
        )
        ws["variance_test_passed"] = variance_preserved

    if show_diagnostic_plot:
        fig, ax = plot_spectra(
            Guu,
            Gvv,
            Guv,
            Gpp,
            frequencies,
            first_frequency_cutoff,
            ff,
            last_frequency_cutoff,
            lf,
            infra_gravity_cutoff,
            ig,
        )
        ws["figure"] = fig
        ws["axis"] = ax

    return ws


def plot_spectra(
    Guu,
    Gvv,
    Guv,
    Gpp,
    frequencies,
    first_frequency_cutoff,
    ff,
    last_frequency_cutoff,
    lf,
    infra_gravity_cutoff,
    ig,
):
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
    ax[0].plot(frequencies, Guv, label="Guv")
    ax[0].plot(frequencies, Guu, label="Guu")
    ax[0].plot(frequencies, Gvv, label="Gvv")
    ylims = ax[0].get_ylim()
    ax[0].plot(
        [frequencies[ff], frequencies[ff]],
        ylims,
        label="ff = %3.1f s" % (1 / frequencies[ff]),
    )
    ax[0].plot(
        [frequencies[lf], frequencies[lf]],
        ylims,
        label="lf = %3.1f s" % (1 / frequencies[lf]),
    )
    ax[0].plot(
        [frequencies[ig], frequencies[ig]],
        ylims,
        label="ig = %3.1f s" % (1 / frequencies[ig]),
    )
    plt.text(
        0.5,
        0.9,
        "first_frequency_cutoff is {} found at #{} = {} Hz".format(
            first_frequency_cutoff, ff, frequencies[ff]
        ),
        transform=ax[0].transAxes,
    )
    plt.text(
        0.5,
        0.8,
        "last_frequency_cutoff is {} found at #{} = {} Hz".format(
            last_frequency_cutoff, lf, frequencies[lf]
        ),
        transform=ax[0].transAxes,
    )
    plt.text(
        0.5,
        0.7,
        "infra-gravity cutoff is {} found at #{} = {} Hz".format(
            infra_gravity_cutoff, ig, frequencies[ig]
        ),
        transform=ax[0].transAxes,
    )
    ax[0].legend()
    ax[1].plot(frequencies, Gpp, label="Gpp")
    ylims = ax[1].get_ylim()
    ax[1].plot(
        [frequencies[ff], frequencies[ff]],
        ylims,
        label="ff = %3.1f s" % (1 / frequencies[ff]),
    )
    ax[1].plot(
        [frequencies[lf], frequencies[lf]],
        ylims,
        label="lf = %3.1f s" % (1 / frequencies[lf]),
    )
    ax[1].plot(
        [frequencies[ig], frequencies[ig]],
        ylims,
        label="ig = %3.1f s" % (1 / frequencies[ig]),
    )
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

    print(
        "These pairs of time-domain and spectral stats should be very close to equal:"
    )
    varp = np.var(p)
    varP = np.sum(Gpp * df)
    print("var(p) {} == sum(Gpp*df) {}".format(varp, varP))
    if np.abs((varp - varP) / varp * 100) > allowable_error:
        result = False
    varu = np.var(u)
    varU = np.sum(Guu * df)
    print("var(u) {} == sum(Guu*df) {}".format(varu, varU))
    if np.abs((varu - varU) / varu * 100) > allowable_error:
        result = False
    varv = np.var(v)
    varV = np.sum(Gvv * df)
    print("var(v) {} == sum(Gvv*df) {}".format(varv, varV))
    if np.abs((varv - varV) / varv * 100) > allowable_error:
        result = False
    print(
        "sqrt(varv) {} == sqrt(mean(v^2) {}".format(
            np.sqrt(varv), np.sqrt(np.mean(v**2))
        )
    )
    varuv = varu + varv
    varUV = varU + varV
    print(
        "sqrt(2.*(varu+varv)) {} == ".format(np.sqrt(2.0 * varuv))
        + "np.sqrt(2*(np.sum(Guu*df) + np.sum(Gvv*df))) {}".format(np.sqrt(2.0 * varUV))
    )
    percent_error = np.abs(
        (np.sqrt(2.0 * varuv) - np.sqrt(2.0 * varUV)) / np.sqrt(2.0 * varuv) * 100
    )
    print(f"percent_error = {percent_error}")
    print(f"allowable_error = {allowable_error}")
    # careful here, it seems that one can't test for 0.0 == 0.0 via the variables
    if (allowable_error == 0.0) & (percent_error == 0.0):
        result = True
    elif (allowable_error > 0.0) & (percent_error < allowable_error):
        result = True

    return result


def puv_qaqc(ds):
    bads = ds["puv_Hrmsu_tail"].isnull()
    for k in ["puv_phir_tail", "puv_azr_tail"]:
        ds[k][bads] = np.nan

    return ds


def next_power_of_2(x):
    return 1 if x == 0 else 2 ** (x - 1).bit_length()


def floor_power_of_2(x):
    return 1 if x == 0 else 2 ** ((x).bit_length() - 1)


def call_puv_quick_vectorized(
    ds,
    pvar=None,
    u=None,
    v=None,
    psh=None,
    huv=None,
    first_frequency_cutoff=None,
    last_frequency_cutoff=None,
):
    desc = {
        "Hrmsp": f"Hrms (=Hmo) from pressure in freq. band {first_frequency_cutoff} <= f <= {last_frequency_cutoff}",
        "Hrmsu": f"Hrms from u,v in freq. band {first_frequency_cutoff} <= f <= {last_frequency_cutoff}",
        "ubr": f"Representative orbital velocity amplitude in freq. band {first_frequency_cutoff} <= f <= {last_frequency_cutoff}",
        "omegar": f"Representative orbital velocity (radian frequency) in freq. band {first_frequency_cutoff} <= f <= {last_frequency_cutoff}",
        "Tr": f"Representative orbital velocity period in freq. band {first_frequency_cutoff} <= f <= {last_frequency_cutoff}",
        "Tpp": f"Peak period from pressure in freq. band {first_frequency_cutoff} <= f <= {last_frequency_cutoff}",
        "Tpu": f"Peak period from velocity in freq. band {first_frequency_cutoff} <= f <= {last_frequency_cutoff}",
        "phir": "Representative orbital velocity direction (angles from x-axis, positive ccw)",
        "azr": "Representative orb. velocity direction (deg; geographic azimuth; ambiguous =/- 180 degrees)",
        "ublo": f"ubr in freq. band f <= {first_frequency_cutoff}",
        "ubhi": f"ubr in freq. band f >= {last_frequency_cutoff}",
        "ubig": f"ubr in infra-gravity freq. band {first_frequency_cutoff} f <= 1/20",
        "Hrmsp_tail": "Hrms (=Hmo) from pressure with f^-4 tail applied",
        "Hrmsu_tail": "Hrms from u,v with f^-4 tail applied",
        "phir_tail": "Representative orbital velocity direction (angles from x-axis, positive ccw) with f^-4 tail applied",
        "azr_tail": "Representative orb. velocity direction (deg; geographic azimuth; ambiguous =/- 180 degrees) with f^-4 tail applied",
        "Snp": f"Pressure-derived non-directional wave energy spectrum in freq. band {first_frequency_cutoff} <= f <= {last_frequency_cutoff}",
        "Snp_tail": "Pressure-derived non-directional wave energy spectrum with f^-4 tail applied",
        "Snu": f"Velocity-derived non-directional wave energy spectrum in freq. band {first_frequency_cutoff} <= f <= {last_frequency_cutoff}",
        "Snu_tail": "Velocity-derived non-directional wave energy spectrum with f^-4 tail applied",
        "frequencies": "Frequency",
        "fclip": "Frequency",
        "Guu": "u-velocity spectrum",
        "Gvv": "v-velocity spectrum",
        "Guv": "uv-velocity spectrum",
        "Gpp": "pressure spectrum",
    }
    standard_name = {
        "Tpp": "sea_surface_wave_period_at_variance_spectral_density_maximum",
        "Tpu": "sea_surface_wave_period_at_variance_spectral_density_maximum",
        "phir": "sea_surface_wave_from_direction",
        "phir_tail": "sea_surface_wave_from_direction",
        "azr": "sea_surface_wave_from_direction",
        "azr_tail": "sea_surface_wave_from_direction",
        "Snp": "sea_surface_wave_variance_spectral_density",
        "Snp_tail": "sea_surface_wave_variance_spectral_density",
        "Snu": "sea_surface_wave_variance_spectral_density",
        "Snu_tail": "sea_surface_wave_variance_spectral_density",
    }
    unit = {
        "Hrmsp": "m",
        "Hrmsu": "m",
        "ubr": "m s-1",
        "omegar": "rad s-1",
        "Tr": "s",
        "Tpp": "s",
        "Tpu": "s",
        "phir": "radians",
        "phir_tail": "radians",
        "azr": "degrees",
        "azr_tail": "degrees",
        "ublo": "m s-1",
        "ubhi": "m s-1",
        "ubig": "m s-1",
        "Hrmsp_tail": "m",
        "Hrmsu_tail": "m",
        "Snp": "m2/Hz",
        "Snp_tail": "m2/Hz",
        "Snu": "m2/Hz",
        "Snu_tail": "m2/Hz",
    }

    # puvs = {k: np.full_like(ds["time"].values, np.nan, dtype=float) for k in desc}
    # puvs = {k: [] for k in desc}

    puvs = puv_quick_vectorized(
        ds[pvar].squeeze(),
        u,
        v,
        ds[pvar].squeeze().mean(dim="sample").values + psh,
        psh,
        huv,
        1 / ds.attrs["sample_interval"],
        first_frequency_cutoff=first_frequency_cutoff,
        last_frequency_cutoff=last_frequency_cutoff,
    )

    ds["puv_frequency"] = xr.DataArray(
        puvs["frequencies"],
        dims="puv_frequency",
        attrs={"standard_name": "sea_surface_wave_frequency", "units": "Hz"},
    )
    ds["puv_frequency_clipped"] = xr.DataArray(
        puvs["fclip"],
        dims="puv_frequency_clipped",
        attrs={"standard_name": "sea_surface_wave_frequency", "units": "Hz"},
    )

    for k in puvs:
        if k == "frequencies" or k == "fclip":
            continue
        if puvs[k].ndim == 1:
            ds["puv_" + k] = xr.DataArray(puvs[k], dims="time")
        elif puvs[k].ndim == 2:
            try:
                ds["puv_" + k] = xr.DataArray(puvs[k], dims=["time", "puv_frequency"])
            except ValueError:
                ds["puv_" + k] = xr.DataArray(
                    puvs[k], dims=["time", "puv_frequency_clipped"]
                )
        if k in desc:
            ds["puv_" + k].attrs["long_name"] = desc[k]
        if k in standard_name:
            ds["puv_" + k].attrs["standard_name"] = standard_name[k]
        if k in unit:
            ds["puv_" + k].attrs["units"] = unit[k]

    ds["puv_Hsp"] = np.sqrt(2) * ds["puv_Hrmsp"]
    ds["puv_Hsp"].attrs["long_name"] = "Hs computed via sqrt(2) * Hrmsp"
    ds["puv_Hsp"].attrs["standard_name"] = "sea_surface_wave_significant_height"
    ds["puv_Hsp"].attrs["units"] = "m"

    ds["puv_Hsu"] = np.sqrt(2) * ds["puv_Hrmsu"]
    ds["puv_Hsu"].attrs["long_name"] = "Hs computed via sqrt(2) * Hrmsu"
    ds["puv_Hsu"].attrs["standard_name"] = "sea_surface_wave_significant_height"
    ds["puv_Hsu"].attrs["units"] = "m"

    ds["puv_Hsp_tail"] = np.sqrt(2) * ds["puv_Hrmsp_tail"]
    ds["puv_Hsp_tail"].attrs[
        "long_name"
    ] = "Hs computed via sqrt(2) * Hrmsp with f^-4 tail applied"
    ds["puv_Hsp_tail"].attrs["standard_name"] = "sea_surface_wave_significant_height"
    ds["puv_Hsp_tail"].attrs["units"] = "m"

    ds["puv_Hsu_tail"] = np.sqrt(2) * ds["puv_Hrmsu_tail"]
    ds["puv_Hsu_tail"].attrs[
        "long_name"
    ] = "Hs computed via sqrt(2) * Hrmsu with f^-4 tail applied"
    ds["puv_Hsu_tail"].attrs["standard_name"] = "sea_surface_wave_significant_height"
    ds["puv_Hsu_tail"].attrs["units"] = "m"

    return ds


def var_wave_burst_fill_nans(ds, var):
    """Fill NaN values in wave burst data if within tolerance and less than 10% of values are NaN"""

    # check for NaNs and fill if possible
    if var.isnull().any():

        nsamps = len(var)

        # Check that number of samples to fill is less than 10% of samples in burst
        if var.isnull().sum().values < nsamps * 0.1:
            # set tolerance for filling gaps in wave burst data
            if "wavedat_tolerance" not in ds.attrs:
                ds.attrs["wavedat_tolerance"] = "2 s"

            # convert tolerance to samples because filling across sample dimension
            tolsamp = (
                int(ds.attrs["wavedat_tolerance"].split()[0]) * ds.attrs["sample_rate"]
            )

            print(
                f" Fill NaNs in {var.name} burst if within tolerance of {tolsamp} samples"
            )

            nnans = var.isnull().sum().values
            var = var.where(~var.isnull().compute(), drop=True)

            var = var.reindex(sample=range(nsamps), method="nearest", tolerance=tolsamp)
            filled_nans = nnans - var.isnull().sum().values
            print(
                f"Filled {filled_nans} of {nnans} NaNs in {var.name} for wave burst at {pd.to_datetime(var.time.values):%Y-%m-%d %H:%M:%S}"
            )

        else:
            print(
                f" Not attempting to fill NaNs ({var.isnull().sum().values}), exceeds 10% of samples ({nsamps})"
            )

    return var


def make_wave_bursts_mi(
    ds,
    wave_vars=[
        "sample",
        "P_1",
        "P_1ac",
        "Tx_1211",
        "brangeAST",
        "ast_quality",
        "u_1205",
        "v_1206",
        "w_1204",
        "vel",
        "cor",
        "amp",
    ],
):
    """
    Reshape CONTINUOUS data set into burst shape using multi-indexing for purpose of wave analysis
    """
    for k in ds.data_vars:
        if k not in wave_vars:
            ds = ds.drop_vars(k)

    # average_interval is [sec] interval for wave statistics for continuous data
    ds.attrs["wave_samples_per_burst"] = int(
        ds.attrs["wave_interval"] / ds.attrs["sample_interval"]
    )
    nsamps = ds.attrs["wave_samples_per_burst"]

    if len(ds.time) % nsamps != 0:
        ds = ds.isel(time=slice(0, -(int(len(ds.time) % (nsamps)))))

    # create samples & new_time
    x = np.arange(nsamps)
    y = ds["time"][0:-1:nsamps]

    # make new arrays for multi-index
    samp, t = np.meshgrid(x, y)
    s = samp.flatten()
    t3 = t.flatten()

    # create multi-index
    ind = pd.MultiIndex.from_arrays((t3, s), names=("new_time", "new_sample"))
    # unstack to make burst shape dataset
    ds = ds.assign_coords(
        xr.Coordinates.from_pandas_multiindex(ind, dim="time")
    ).unstack("time")
    # dsa = dsa.unstack()
    if "sample" in ds.data_vars:
        ds = ds.drop_vars("sample").rename({"new_time": "time", "new_sample": "sample"})
    else:
        ds = ds.rename({"new_time": "time", "new_sample": "sample"})

    return ds


def clip_f(spec, nsamps=None, fs=None):
    """Clip frequencies in wave spectra data set"""

    if nsamps and fs:
        nyfreq = float(fs) / 2
        flo = np.round(
            1 / (nsamps / fs / 32.0), 3
        )  # set minimum frequency to length of wave burst samples use divided by 32
        if nyfreq > 2:
            fhi = 2
        else:
            fhi = nyfreq

        # apply buffer to each end of slice to include end points
        spec = spec.sel(frequency=slice(flo - 0.005, fhi + 0.005))
    else:
        print(
            "Cannot clip wave spectral frequency because inputs nsamps and/or fs are missing"
        )

    return spec
