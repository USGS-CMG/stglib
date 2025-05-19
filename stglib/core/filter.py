import numpy as np
import scipy.signal as spsig
import xarray as xr

from . import utils


def butter_filt(sig, sr, cutfreq, ftype, ford=4):
    """
    butterworth filter using sosfiltfilt in scipy.signal

    Parameters
    ----------
        sig - signal to be filtered (array)
        sr = sample rate of signal (Hz)
        cutfreq = cutoff frequency for filter (Hz) length = 1 or 2 (for bandpass)
        ftype = type of filter, options = ['lowpass', 'highpass', 'bandpass']
        ford = filter order (default = 4)

    Returns
    -------
        filtered signal using specified order (default = 4) butterworth filter

    """
    sos = spsig.butter(ford, cutfreq, btype=ftype, fs=sr, output="sos")

    return spsig.sosfiltfilt(sos, sig)


def make_butter_filt(ds, var, sr, cutfreq, ftype):
    """
    Create smoothed data using specified butterworth filter type, order, and cutoff for user specified variables

    Parameters
    ----------
        ds - xarray dataset that contains user specified variable
        var - user specified variable
        sr - sample rate of data (hertz)
        cutfreq - cutoff frequency(s) for filter (hertz)
        ftype - user specified filter type (lowpass, highpass, or bandpass)

    Returns
    -------
        ds - dataset with specified variable smoothed/filtered with the specified butterworth filter type, order, and cutoff

    """
    if "filter_order" in ds.attrs:
        ford = ds.attrs["filter_order"]
    else:
        ford = 4

    if ds[var].ndim == 1 and "time" in ds[var].dims:
        print(f"Applying {ftype} filter to {var}")
        filtered = butter_filt(ds[var].values, sr, cutfreq, ftype, ford)
        ds[var][:] = filtered

        notetxt = f"Values filtered using order = {ford} butterworth {ftype} filter with {cutfreq} cutoff frequency. "
        ds = utils.insert_note(ds, var, notetxt)

    elif ds[var].ndim == 2 and "time" in ds[var].dims and "sample" in ds[var].dims:
        print(f"Applying {ftype} filter to {var} burst data")
        for k in ds["time"]:

            filtered = butter_filt(ds[var].sel(time=k).values, sr, cutfreq, ftype, ford)
            ds[var].loc[dict(time=k)] = filtered

        notetxt = f"Values filtered using order = {ford} butterworth {ftype} filter with {cutfreq} Hz cutoff frequency. "
        ds = utils.insert_note(ds, var, notetxt)

    elif ds[var].ndim == 2 and "time" in ds[var].dims and "z" in ds[var].dims:
        print(f"Applying {ftype} filter to {var} profile data")
        for k in ds["z"]:

            filtered = butter_filt(ds[var].sel(z=k).values, sr, cutfreq, ftype, ford)
            ds[var].loc[dict(z=k)] = filtered

        notetxt = f"Values filtered using order = {ford} butterworth {ftype} filter with {cutfreq} Hz cutoff frequency. "
        ds = utils.insert_note(ds, var, notetxt)

    elif ds[var].ndim == 2 and "time" in ds[var].dims and "bindist" in ds[var].dims:
        print(f"Applying {ftype} filter to {var} profile data")
        for k in ds["bindist"]:

            filtered = butter_filt(
                ds[var].sel(bindist=k).values, sr, cutfreq, ftype, ford
            )
            ds[var].loc[dict(bindist=k)] = filtered

        notetxt = f"Values filtered using order = {ford} butterworth {ftype} filter with {cutfreq} Hz cutoff frequency. "
        ds = utils.insert_note(ds, var, notetxt)

    else:
        print(
            f"Not able to apply {ftype} filter because only 'time', ('time','sample'), ('time','z'), ('time', 'bindist') dimensions are handled and {var} dims are {ds[var].dims}"
        )

    return ds


def apply_butter_filt(ds, var):
    """
    Construct and call butterworth filter from user specified config.yaml file

    Parameters
    ----------
        ds - xarray dataset with user specified variable
        var - user specified variable

    Returns
    -------
        ds - dataset with specified variable smoothed/filtered with the specified butterworth filter type, order, and cutoff

    """
    if (
        var + "_lowpass_filt" in ds.attrs
        or var + "_highpass_filt" in ds.attrs
        or var + "_bandpass_filt" in ds.attrs
    ):
        ds[var + "_unfiltered"] = ds[var].copy()
        if "long_name" in ds[var + "_unfiltered"].attrs:
            ds[var + "_unfiltered"].attrs["long_name"] = (
                ds[var + "_unfiltered"].attrs["long_name"] + " (unfiltered)"
            )

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
                f"sample_rate or sample _interval do not exist in global attributes, cannot apply lowpass filter to {var}. "
            )

    return ds


def apply_med_filt(ds, var):
    """
    Construct and apply N-point median filter to user specified variable and N

    Parameters
    ----------
        ds - xarray dataset with user specified variable
        var - user specified variable

    Returns
    -------
        ds - dataset with user specified variable smoothed/filtered with the user specified N points (kernel size).
    """
    if var + "_med_filt" in ds.attrs:

        if var + "_unfiltered" not in ds.data_vars:
            ds[var + "_unfiltered"] = ds[var].copy()
            if "long_name" in ds[var + "_unfiltered"].attrs:
                ds[var + "_unfiltered"].attrs["long_name"] = (
                    ds[var + "_unfiltered"].attrs["long_name"] + " (unfiltered)"
                )

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

            elif ds[var].ndim == 2 and "time" in ds[var].dims and "z" in ds[var].dims:
                print(
                    f"Applying {kernel_size} point median filter to {var} profile data"
                )

                for k in ds["z"]:

                    filtered = spsig.medfilt(ds[var].sel(z=k).values, kernel_size)
                    ds[var].loc[dict(z=k)] = filtered

                notetxt = f"Values filtered using {kernel_size} point median filter. "
                ds = utils.insert_note(ds, var, notetxt)

            else:
                print(
                    f"Not able to apply median filter because only 'time', ('time','sample'), or ('time', 'z') dimensions are handled and {var} dims are {ds[var].dims}"
                )

        else:
            raise ValueError(
                f"Not able to apply median filter because kernel size specified {kernel_size} is not an odd whole number"
            )

    return ds


def filter_vel(ds, data_vars=["U", "V", "W"]):
    """Filter velocity variables if specified using Butterworth filter - save orignal unfiltered data also"""

    if "vel_filter_method" in ds.attrs and "vel_filter_cut" in ds.attrs:
        # loop through data variables
        if "sample_rate" in ds.attrs or "sample_interval" in ds.attrs:
            # check to make sure sample_rate or sample_intreval attributes exits.
            if "sample_rate" in ds.attrs:
                sr = ds.attrs["sample_rate"]
            else:
                sr = 1 / ds.attrs["sample_interval"]

            # loop through data variables
            for var in data_vars:

                if var in ds.data_vars:

                    # back-up attrs
                    attrsbak = ds[var].attrs
                    ds[var + "_unfiltered"] = ds[var].copy()
                    ds[var + "_unfiltered"].attrs = attrsbak

                    ftype = ds.attrs["vel_filter_method"]

                    if ftype == "lowpass" or ftype == "highpass":
                        cutfreq = 1 / ds.attrs["vel_filter_cut"]
                        ds = make_butter_filt(ds, var, sr, cutfreq, ftype)

                    elif ftype == "bandpass":
                        cutfreq_lo = 1 / ds.attrs["vel_filter_cut"][0]
                        cutfreq_hi = 1 / ds.attrs["vel_filter_cut"][1]
                        cutfreq = [cutfreq_lo, cutfreq_hi]
                        ds = make_butter_filt(ds, var, sr, cutfreq, ftype)
                    else:
                        raise ValueError(
                            f"filter method specified {ftype} is not in avaiable list [lowpass, highpass, or bandpass]."
                        )

        else:
            raise ValueError(
                f"sample_rate or sample _interval do not exist in global attributes, cannot apply {ds.attrs['vel_filter_method']} filter to velocity variables. "
            )

    return ds
