import numpy as np
import scipy.signal as spsig
import xarray as xr

from . import utils


def butter_filt(sig, sr, cutfreq, ftype):
    # Returns filtered signal using 2nd order butterworth filter
    # Inputs:
    # sig - signal to be filtered (array)
    # sr = sample rate of signal (Hz)
    # cutfreq = cutoff frequency for filter (Hz) length = 1 or 2 (for bandpass)
    # ftype = type of filter, options = ['lowpass', 'highpass', 'bandpass']

    sos = spsig.butter(2, cutfreq, btype=ftype, fs=sr, output="sos")

    return spsig.sosfiltfilt(sos, sig)


def make_butter_filt(ds, var, sr, cutfreq, ftype):

    # create smoothed data using specified butterworth filter type and cutoff
    if ds[var].ndim == 1 and "time" in ds[var].dims:
        print(f"Applying {ftype} filter to {var}")
        filtered = butter_filt(ds[var].values, sr, cutfreq, ftype)
        ds[var][:] = filtered

        notetxt = f"Values filtered using 2nd order butterworth {ftype} filter with {cutfreq} cutoff frequency. "
        ds = utils.insert_note(ds, var, notetxt)

    elif ds[var].ndim == 2 and "time" in ds[var].dims and "sample" in ds[var].dims:
        print(f"Applying {ftype} filter to {var} burst data")
        for k in ds["time"]:

            filtered = butter_filt(ds[var].sel(time=k).values, sr, cutfreq, ftype)
            ds[var].loc[dict(time=k)] = filtered

        notetxt = f"Values filtered using 2nd order butterworth {ftype} filter with {cutfreq} Hz cutoff frequency. "
        ds = utils.insert_note(ds, var, notetxt)

    elif ds[var].ndim == 2 and "time" in ds[var].dims and "z" in ds[var].dims:
        print(f"Applying {ftype} filter to {var} profile data")
        for k in ds["z"]:

            filtered = butter_filt(ds[var].sel(z=k).values, sr, cutfreq, ftype)
            ds[var].loc[dict(z=k)] = filtered

        notetxt = f"Values filtered using 2nd order butterworth {ftype} filter with {cutfreq} Hz cutoff frequency. "
        ds = utils.insert_note(ds, var, notetxt)

    else:
        print(
            f"Not able to apply {ftype} filter because only 'time' , ('time','sample'), or ('time','z') dimensions are handled and {var} dims are {ds[var].dims}"
        )

    return ds


def apply_butter_filt(ds, var):

    if (
        var + "_lowpass_filt" in ds.attrs
        or var + "_highpass_filt" in ds.attrs
        or var + "_bandpass_filt" in ds.attrs
    ):

        if (
            "sample_rate" in ds.attrs or "sample_interval" in ds.attrs
        ):  # check to make sure sample_rate or sample_intreval attributes exits.
            if "sample_rate" in ds.attrs:
                sr = ds.attrs["sample_rate"]
            else:
                sr = 1 / ds.attrs["sample_interval"]

            if var + "_lowpass_filt" in ds.attrs:
                ftype = "lowpass"
                cutfreq = 1 / ds.attrs[var + "_lowpass_filt"]
                ds = make_butter_filt(ds, var, sr, cutfreq, ftype)

            elif var + "_highpass_filt" in ds.attrs:
                ftype = "highpass"
                cutfreq = 1 / ds.attrs[var + "_highpass_filt"]
                ds = make_butter_filt(ds, var, sr, cutfreq, ftype)

            elif var + "_bandpass_filt" in ds.attrs:
                ftype = "bandpass"
                cutfreq_lo = 1 / ds.attrs[var + "_bandpass_filt"][0]
                cutfreq_hi = 1 / ds.attrs[var + "_bandpass_filt"][1]
                cutfreq = [cutfreq_lo, cutfreq_hi]
                print(cutfreq)
                ds = make_butter_filt(ds, var, sr, cutfreq, ftype)

        else:
            raise ValueError(
                f"sample_rate or sample _interval do not exits in global attributes, can not apply lowpass filter to {var}. "
            )

    return ds


def apply_med_filt(ds, var):
    # print(var)
    if var + "_med_filt" in ds.attrs:

        kernel_size = ds.attrs[var + "_med_filt"]
        # make sure kernel_size is odd number
        if kernel_size % 2 == 1:

            if ds[var].ndim == 1 and "time" in ds[var].dims:
                print(f"Applying {kernel_size} point median filter to {var}")
                filtered = spsig.medfilt(ds[var].values, kernel_size)
                ds[var][:] = filtered

                notetxt = f"Values filtered using {kernel_size} point median filter. "
                ds = utils.insert_note(ds, var, notetxt)

            elif (
                ds[var].ndim == 2
                and "time" in ds[var].dims
                and "sample" in ds[var].dims
            ):
                print(f"Applying {kernel_size} point median filter to {var} burst data")

                for k in ds["time"]:
                    filtered = spsig.medfilt(ds[var].sel(time=k).values, kernel_size)
                    ds[var].loc[dict(time=k)] = filtered

                notetxt = f"Values filtered using {kernel_size} point median filter. "
                ds = utils.insert_note(ds, var, notetxt)

            else:
                print(
                    f"Not able to apply median filter because only 'time' or ('time','sample') dimensions are handled and {var} dims are {ds[var].dims}"
                )

        else:
            raise ValueError(
                f"Not able to apply median filter because kernel size specified {kernel_size} is not an odd whole number"
            )

    return ds
