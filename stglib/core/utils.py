import csv
import datetime
import inspect
import os
import platform
import re
import sqlite3
import sys
import warnings
from pathlib import Path

import netCDF4
import numpy as np
import pandas as pd
import scipy.io as spio
import xarray as xr

import stglib

from . import filter


def is_cf(ds):
    if ("Conventions" in ds.attrs) and ("CF" in ds.attrs["Conventions"]):
        return True
    else:
        return False


def ensure_cf(ds):
    if not is_cf(ds):
        raise ValueError(
            "Non-CF Conventions are not supported. Ensure you are setting Conventions appropriately."
        )

    if ds.attrs["Conventions"] != "CF-1.9":
        warnings.warn(
            f"You are using a version of the CF Conventions ({ds.attrs['Conventions']}) that is not the latest supported version (CF-1.9). Consider changing to CF-1.9."
        )

    return ds


def check_compliance(nc_file, conventions="CF-1.9"):
    from compliance_checker.runner import CheckSuite, ComplianceChecker

    check_suite = CheckSuite()
    check_suite.load_all_available_checkers()

    verbose = 1
    criteria = "normal"
    output_filename = f"{nc_file}.cfcheck.txt"
    output_format = "text"
    checker_names = [conventions.lower().replace("-", ":")]

    print(
        f"*** Checking CF compliance. Displaying results of compliance checker below (also saved in {output_filename})."
    )

    return_value, errors = ComplianceChecker.run_checker(
        str(nc_file),
        checker_names,
        verbose,
        criteria,
        output_filename=output_filename,
        output_format=output_format,
    )

    with open(output_filename) as f:
        print(f.read())


def clip_ds(ds, wvs=False):
    """
    Clip an xarray Dataset from metadata, either via good_ens or
    Deployment_date and Recovery_date.

    wvs only applies to Aquadopp waves here. It is distinct from waves flag
    because AQD waves can have a different sampling interval than AQD currents
    """

    print(
        "first burst in full file: {}, idx {}".format(
            ds["time"].min().values, np.argmin(ds["time"].values)
        )
    )
    print(
        "last burst in full file: {}, idx {}".format(
            ds["time"].max().values, np.argmax(ds["time"].values)
        )
    )

    # clip either by ensemble indices or by the deployment and recovery
    # date specified in metadata
    if "good_ens" in ds.attrs and not wvs:
        # we have good ensemble indices in the metadata

        # so we can deal with multiple good_ens ranges, or just a single range
        good_ens = ds.attrs["good_ens"]
        goods = []

        for n in range(0, len(good_ens), 2):
            goods.append(np.arange(good_ens[n], good_ens[n + 1]))
        goods = np.hstack(goods)
        ds = ds.isel(time=goods)

        histtext = "Data clipped using good_ens values of {}.".format(str(good_ens))

        ds = insert_history(ds, histtext)

    elif "good_ens_wvs" in ds.attrs and wvs:
        good_ens = ds.attrs["good_ens_wvs"]

        goods = np.arange(good_ens[0], good_ens[1])

        ds = ds.isel(time=goods)

        histtext = "Data clipped using good_ens_wvs values of {}.".format(str(good_ens))

        ds = insert_history(ds, histtext)

    elif "good_dates" in ds.attrs:
        # clip by start/end dates that are not Deployment_date and Recovery_date

        good_dates = ds.attrs["good_dates"]
        goods = []

        for n in range(0, len(good_dates), 2):
            toappend = (ds.time > np.datetime64(good_dates[n])) & (
                ds.time <= np.datetime64(good_dates[n + 1])
            )
            where = np.where(toappend)[0]
            print(f"{good_dates[n]=}, {where.min()=}")
            print(f"{good_dates[n+1]=}, {where.max()=}")
            goods.append(toappend)

        goods = np.any(goods, axis=0)
        ds = ds.isel(time=goods)

        # where = np.where(
        #     (ds["time"].values >= np.datetime64(ds.attrs["good_dates"][0]))
        #     & (ds["time"].values <= np.datetime64(ds.attrs["good_dates"][1]))
        # )[0]
        # print("good_dates[0] {}, idx {}".format(ds.attrs["good_dates"][0], where.min()))
        # print("good_dates[1] {}, idx {}".format(ds.attrs["good_dates"][1], where.max()))
        # ds = ds.sel(time=slice(ds.attrs["good_dates"][0], ds.attrs["good_dates"][1]))

        histtext = "Data clipped using good_dates of {}.".format(
            ds.attrs["good_dates"],
        )

        ds = insert_history(ds, histtext)

    elif "Deployment_date" in ds.attrs and "Recovery_date" in ds.attrs:
        # we clip by the times in/out of water as specified in the metadata

        ds = ds.sel(time=slice(ds.attrs["Deployment_date"], ds.attrs["Recovery_date"]))

        histtext = (
            "Data clipped using Deployment_date of {} and Recovery_date of {}."
        ).format(
            ds.attrs["Deployment_date"],
            ds.attrs["Recovery_date"],
        )

        ds = insert_history(ds, histtext)
    else:
        # do nothing
        print("Did not clip data; no values specified in metadata")

    try:
        print("first burst in trimmed file:", ds["time"].min().values)
        print("last burst in trimmed file:", ds["time"].max().values)
    except ValueError:
        raise (
            ValueError(
                "No valid time values in trimmed dataset. Are you "
                "sure you sure you specified Deployment and Recovery "
                "dates correctly?"
            )
        )

    return ds


def add_min_max(ds, exclude_vars=None):
    """
    Add minimum and maximum values to variables in NC or CDF files
    This function assumes the data are in xarray DataArrays within Datasets
    """

    exclude = list(ds.dims)
    [exclude.append(k) for k in ds.variables if re.match("time*", k)]
    exclude.extend(["TIM", "TransMatrix", "orientmat"])
    if exclude_vars:
        exclude.extend(exclude_vars)

    alloweddims = [
        "time",
        "sample",
        "depth",
        "z",
        "x",
        "y",
        "sweep",
        "scan",
        "points",
        "frequency",
        "bin_along",
        "bin_across",
        "direction",
        "diwasp_frequency",
        "diwasp_direction",
        "puv_frequency",
        "puv_frequency_clipped",
        "beam",
        "inst",
        "earth",
        "inst4",
        "earth4",
        "q",
        "depthsen",
        "zsen",
    ]

    for k in ds.variables:
        if k not in exclude:
            kwargs = {"dim": tuple(d for d in alloweddims if d in ds[k].dims)}

            ds[k].attrs.update(
                {
                    "minimum": ds[k].min(**kwargs).squeeze().values,
                    "maximum": ds[k].max(**kwargs).squeeze().values,
                }
            )

    return ds


def insert_history(ds, histtext):
    toinsert = "{}: {}\n".format(
        datetime.datetime.now(datetime.timezone.utc).isoformat(), histtext
    )

    print(toinsert.rstrip())

    if "history" in ds.attrs:
        ds.attrs["history"] = ds.attrs["history"] + toinsert
    else:
        ds.attrs["history"] = toinsert

    return ds


def add_history(ds):
    if is_cf(ds):
        histtext = "Processed to {} using {}.".format(
            ds.attrs["Conventions"],
            os.path.basename(sys.argv[0]),
        )
    else:
        histtext = "Processed to EPIC using {}.".format(os.path.basename(sys.argv[0]))

    return insert_history(ds, histtext)


def ds_add_waves_history(ds):
    histtext = (
        "Wave statistics computed using scipy.signal.welch(), "
        "assigning cutoff following Jones & Monismith (2007), and "
        "applying f^-4 tail past cutoff."
    )

    return insert_history(ds, histtext)


def ds_add_diwasp_history(ds):
    """
    Add history indicating DIWASP has been applied
    """

    histtext = "Wave statistics computed using DIWASP 1.1GD."

    return insert_history(ds, histtext)


def ds_add_pydiwasp_history(ds):
    """
    Add history indicating DIWASP has been applied
    """

    histtext = f"Directional Wave statistics computed using pyDIWASP with {ds.attrs['diwasp']} input data"

    return insert_history(ds, histtext)


def ds_coord_no_fillvalue(ds):
    for var in [
        "latitude",
        "longitude",
        "depth",
        "time",
        "time2",
        "time_cf",
        "time_2d",
        "time_cf_2d",
        "epic_time",
        "epic_time2",
        "epic_time_2d",
        "epic_time2_2d",
        "TransMatrix",
        "direction",
        "sample",
        "frequency",
        "z",
        "x",
        "y",
        "bindist",
        "zsen",
        "depthsen",
        "sweep",
        "scan",
        "points",
    ]:
        if var in ds:
            ds[var].encoding["_FillValue"] = None

    return ds


def add_standard_names(ds):
    standard_names = {
        "Hdg_1215": "platform_orientation",
        "Ptch_1216": "platform_pitch",
        "Roll_1217": "platform_roll",
        "P_1": "sea_water_pressure",
        # "Tx_1211": "sea_water_temperature",
        "P_1ac": "sea_water_pressure_due_to_sea_water",
        "u_1205": "eastward_sea_water_velocity",
        "v_1206": "northward_sea_water_velocity",
        "w_1204": "upward_sea_water_velocity",
        "SV_80": "speed_of_sound_in_sea_water",
        "AGC_1202": "signal_intensity_from_multibeam_acoustic_doppler_velocity_sensor_in_sea_water",
        "cor_avg": "beam_consistency_indicator_from_multibeam_acoustic_doppler_velocity_profiler_in_sea_water",
    }

    for k in standard_names:
        if k in ds:
            # if "standard_name" not in ds[k].attrs:
            ds[k] = check_update_attrs(ds[k], "standard_name", standard_names[k])
            # ds[k].attrs["standard_name"] = standard_names[k]

    if "Tx_1211" in ds:
        if ds["Tx_1211"].attrs["units"] == "C":
            ds["Tx_1211"].attrs["units"] = "degree_C"

    for n in [1, 2]:
        if f"AnalogInput{n}" in ds:
            for v in [
                "standard_name",
                "long_name",
                "units",
                "comment",
                "institution",
                "references",
                "source",
            ]:
                if f"AnalogInput{n}_{v}" in ds.attrs:
                    ds[f"AnalogInput{n}"] = check_update_attrs(
                        ds[f"AnalogInput{n}"], v, ds.attrs[f"AnalogInput{n}_{v}"]
                    )

    # ds["feature_type_instance"] = xr.DataArray(
    #    [f"{ds.attrs['MOORING']}aqd"], dims="feature_type_instance"
    # )
    # ds["feature_type_instance"].attrs["cf_role"] = "timeseries_id"

    # for k in ds.data_vars:
    #    ds[k].attrs["coverage_content_type"] = "physicalMeasurement"

    return ds


def ds_add_wave_attrs(ds):
    """
    Add wave attributes to variables
    """

    # Update attributes for EPIC and STG compliance
    ds = ds_coord_no_fillvalue(ds)

    ds["time"].attrs.update(
        {"standard_name": "time", "axis": "T", "long_name": "time (UTC)"}
    )

    ds["time"].encoding["dtype"] = "i4"

    if "burst" in ds:
        check_fits_in_int32(ds, "burst")
        ds["burst"].encoding["dtype"] = "i4"
        ds["burst"].attrs["units"] = "1"
        ds["burst"].attrs["long_name"] = "Burst number"

    def add_attributes(var, dsattrs):
        var.attrs.update(
            {
                # "serial_number": dsattrs["serial_number"],
                # "initial_instrument_height": dsattrs["initial_instrument_height"],
                # "nominal_instrument_depth": dsattrs["nominal_instrument_depth"],
                # "height_depth_units": "m",
            }
        )
        # if "INST_TYPE" in dsattrs:
        #    var.attrs["sensor_type"] = dsattrs["INST_TYPE"]

    if "wp_peak" in ds.data_vars:
        ds["wp_peak"].attrs.update(
            {
                "long_name": "Dominant (peak) wave period",
                "units": "s",
                "epic_code": 4063,
                "standard_name": "sea_surface_wave_period_at_variance_spectral_density_maximum",
            }
        )

    if "wp_4060" in ds.data_vars:
        ds["wp_4060"].attrs.update(
            {
                "long_name": "Average wave period",
                "units": "s",
                "epic_code": 4060,
                "standard_name": "sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment",
            }
        )

    if "wh_4061" in ds.data_vars:
        ds["wh_4061"].attrs.update(
            {
                "long_name": "Significant wave height",
                "units": "m",
                "epic_code": 4061,
                "standard_name": "sea_surface_wave_significant_height",
            }
        )

    if "pspec" in ds.data_vars:
        ds["pspec"].attrs.update(
            {
                "long_name": "Pressure-derived non-directional wave energy spectrum",
                "units": "m^2/Hz",
                "note": "Use caution: all spectra are provisional",
                "standard_name": "sea_surface_wave_variance_spectral_density",
            }
        )

    if "sspec" in ds.data_vars:
        ds["sspec"].attrs.update(
            {
                "long_name": "Surface track derived non-directional wave energy spectrum",
                "units": "m^2/Hz",
                "note": "Use caution: all spectra are provisional",
                "standard_name": "sea_surface_wave_variance_spectral_density",
            }
        )

    if "vspec" in ds.data_vars:
        ds["vspec"].attrs.update(
            {
                "long_name": "Velocity-derived non-directional wave energy spectrum",
                "units": "m^2/Hz",
                "note": "Use caution: all spectra are provisional",
                "standard_name": "sea_surface_wave_variance_spectral_density",
            }
        )

    if "frequency" in ds.coords:
        ds["frequency"].attrs.update(
            {
                "standard_name": "sea_surface_wave_frequency",
                "long_name": "Frequency",
                "units": "Hz",
            }
        )

    if "direction" in ds.coords:
        ds["direction"].attrs.update(
            {
                "long_name": "Direction (from, relative to true north)",
                "units": "degrees",
            }
        )

    if "dspec" in ds.data_vars:
        ds["dspec"].attrs.update(
            {
                "long_name": "Directional wave energy spectrum",
                "units": "m^2/Hz/degree",
                "note": "Use caution: all spectra are provisional",
                "standard_name": "sea_surface_wave_directional_variance_spectral_density",
            }
        )

    if "wvdir" in ds.data_vars:
        ds["wvdir"].attrs.update(
            {
                "long_name": (
                    "Direction of peak period " "(from, relative to true north)"
                ),
                "units": "degrees",
                "note": (
                    "Compass direction from which waves are propagating as "
                    "defined by the direction with the greatest energy at "
                    "the peak period"
                ),
                "standard_name": "sea_surface_wave_from_direction_at_variance_spectral_density_maximum",
            }
        )

    if "dwvdir" in ds.data_vars:
        ds["dwvdir"].attrs.update(
            {
                "long_name": (
                    "Dominant wave direction " "(from, relative to true north)"
                ),
                "units": "degrees",
                "note": (
                    "Compass direction from which waves are propagating as "
                    "defined by the direction band with greatest total "
                    "energy summed over all frequencies"
                ),
            }
        )

    if "wd_4062" in ds.data_vars:
        ds["wd_4062"].attrs.update(
            {
                "long_name": "Mean wave direction",
                "units": "degrees",
                "epic_code": 4062,
                "note": "Compass direction from which waves are propagating",
                "standard_name": "sea_surface_wave_from_direction",
            }
        )

    if "diwasp_frequency" in ds.coords:
        ds["diwasp_frequency"].attrs.update(
            {
                "standard_name": "sea_surface_wave_frequency",
                "long_name": "Frequency",
                "units": "Hz",
            }
        )

    if "diwasp_direction" in ds.coords:
        ds["diwasp_direction"].attrs.update(
            {
                "long_name": "Direction (from, relative to true north)",
                "units": "degrees",
            }
        )

    if "diwasp_tp" in ds.data_vars:
        ds["diwasp_tp"].attrs.update(
            {
                "long_name": "Dominant (peak) wave period from pyDIWASP",
                "units": "s",
                "standard_name": "sea_surface_wave_period_at_variance_spectral_density_maximum",
            }
        )

    if "diwasp_tm" in ds.data_vars:
        ds["diwasp_tm"].attrs.update(
            {
                "long_name": "Average wave period",
                "units": "s",
                "epic_code": 4060,
                "standard_name": "sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment",
            }
        )

    if "diwasp_hs" in ds.data_vars:
        ds["diwasp_hs"].attrs.update(
            {
                "long_name": "Significant wave height from pyDIWASP",
                "units": "m",
                "standard_name": "sea_surface_wave_significant_height",
            }
        )

    if "diwasp_dtp" in ds.data_vars:
        ds["diwasp_dtp"].attrs.update(
            {
                "long_name": (
                    "Direction of peak period "
                    "(from, relative to true north) from pyDIWASP"
                ),
                "units": "degrees",
                "note": (
                    "Compass direction from which waves are propagating as "
                    "defined by the direction with the greatest energy at "
                    "the peak period"
                ),
                "standard_name": "sea_surface_wave_from_direction_at_variance_spectral_density_maximum",
            }
        )

    if "diwasp_dp" in ds.data_vars:
        ds["diwasp_dp"].attrs.update(
            {
                "long_name": (
                    "Dominant wave direction "
                    "(from, relative to true north) from pyDIWASP"
                ),
                "units": "degrees",
                "note": (
                    "Compass direction from which waves are propagating as "
                    "defined by the direction band with greatest total "
                    "energy summed over all frequencies"
                ),
            }
        )

    if "diwasp_dm" in ds.data_vars:
        ds["diwasp_dm"].attrs.update(
            {
                "long_name": "Mean wave direction from pyDIWASP",
                "units": "degrees",
                "note": "Compass direction from which waves are propagating",
                "standard_name": "sea_surface_wave_from_direction",
            }
        )

    if "diwasp_dspec" in ds.data_vars:
        ds["diwasp_dspec"].attrs.update(
            {
                "long_name": "Directional wave energy spectrum from pyDIWASP",
                "units": "m^2/Hz/degree",
                "note": "Use caution: all spectra are provisional",
                "standard_name": "sea_surface_wave_directional_variance_spectral_density",
            }
        )

    if "diwasp_fspec" in ds.data_vars:
        ds["diwasp_fspec"].attrs.update(
            {
                "long_name": "Frequency (non-directional) wave energy spectrum from pyDIWASP",
                "units": "m^2/Hz",
                "note": "Use caution: all spectra are provisional",
                "standard_name": "sea_surface_wave_variance_spectral_density",
            }
        )

    if "diwasp_astspec" in ds.data_vars:
        ds["diwasp_astspec"].attrs.update(
            {
                "long_name": "Acoustic Surface Tracking derived frequency (non-directional) wave energy spectrum from pyDIWASP",
                "units": "m^2/Hz",
                "note": "Use caution: all spectra are provisional",
                "standard_name": "sea_surface_wave_variance_spectral_density",
            }
        )

    if "diwasp_pspec" in ds.data_vars:
        ds["diwasp_pspec"].attrs.update(
            {
                "long_name": "Pressure derived frequency (non-directional) wave energy spectrum from pyDIWASP",
                "units": "m^2/Hz",
                "note": "Use caution: all spectra are provisional",
                "standard_name": "sea_surface_wave_variance_spectral_density",
            }
        )

    if "diwasp_vspec" in ds.data_vars:
        ds["diwasp_vspec"].attrs.update(
            {
                "long_name": "Velocity derived frequency (non-directional) wave energy spectrum from pyDIWASP",
                "units": "m^2/Hz",
                "note": "Use caution: all spectra are provisional",
                "standard_name": "sea_surface_wave_variance_spectral_density",
            }
        )

    for var in [
        "wp_peak",
        "wh_4061",
        "wp_4060",
        "wd_4062",
        "pspec",
        "water_depth",
        "dspec",
        "wvdir",
        "dwvdir",
        "diwasp_hs",
        "diwasp_tp",
        "diwasp_dtp",
        "diwasp_dp",
        "diwasp_dm",
        "diwasp_dspec",
        "diwasp_fspec",
        "diwasp_astspec",
        "diwasp_pspec",
        "diwasp_vspec",
    ]:
        if var in ds.variables:
            add_attributes(ds[var], ds.attrs)
            ds[var].attrs.update(
                {"minimum": ds[var].min().values, "maximum": ds[var].max().values}
            )

    return ds


def trim_max_wp(ds):
    """
    QA/QC
    Trim wave data based on maximum wave period as specified in metadata
    """

    if "wp_max" in ds.attrs:
        for var in ["wp_peak", "wp_4060"]:
            ds[var] = ds[var].where(
                (ds["wp_peak"] < ds.attrs["wp_max"])
                & (ds["wp_4060"] < ds.attrs["wp_max"])
            )

        for var in ["wp_peak", "wp_4060"]:
            notetxt = "Values filled where wp_peak, wp_4060 >= {}.".format(
                ds.attrs["wp_max"]
            )

            ds = insert_note(ds, var, notetxt)

    return ds


def trim_min_wh(ds):
    """
    QA/QC
    Trim wave data based on minimum wave height as specified in metadata
    """

    if "wh_min" in ds.attrs:
        for var in ["wp_peak", "wh_4061", "wp_4060"]:
            ds[var] = ds[var].where(ds["wh_4061"] > ds.attrs["wh_min"])

            notetxt = "Values filled where wh_4061 <= {}.".format(ds.attrs["wh_min"])

            ds = insert_note(ds, var, notetxt)

    return ds


def trim_max_wh(ds):
    """
    QA/QC
    Trim wave data based on maximum wave height as specified in metadata
    """

    if "wh_max" in ds.attrs:
        for var in ["wp_peak", "wh_4061", "wp_4060"]:
            ds[var] = ds[var].where(ds["wh_4061"] < ds.attrs["wh_max"])

            notetxt = "Values filled where wh_4061 >= {}.".format(ds.attrs["wh_max"])

            ds = insert_note(ds, var, notetxt)

    return ds


def trim_wp_ratio(ds):
    """
    QA/QC
    Trim wave data based on maximum ratio of wp_peak to wp_4060
    """

    if "wp_ratio" in ds.attrs:
        for var in ["wp_peak", "wp_4060"]:
            ds[var] = ds[var].where(
                ds["wp_peak"] / ds["wp_4060"] < ds.attrs["wp_ratio"]
            )

        for var in ["wp_peak", "wp_4060"]:
            notetxt = "Values filled where wp_peak:wp_4060 >= {}.".format(
                ds.attrs["wp_ratio"]
            )

            ds = insert_note(ds, var, notetxt)

    return ds


def write_metadata(ds, metadata):
    """Write out all metadata to CDF file"""

    for k in metadata:
        # recursively write out instmeta
        if k == "instmeta":
            for x in metadata[k]:
                ds = check_update_attrs(ds, x, metadata[k][x])
        else:
            ds = check_update_attrs(ds, k, metadata[k])

    for v in ["Deployment_date", "Recovery_date"]:
        if v in ds.attrs:
            # look for errant quotes in dates
            ds.attrs[v] = ds.attrs[v].replace('"', "")

    f = os.path.basename(inspect.stack()[1][1])

    histtext = (
        "Processed using {} with stglib {}, xarray {}, NumPy {}, netCDF4 {}, Python {}."
    ).format(
        f,
        stglib.__version__,
        xr.__version__,
        np.__version__,
        netCDF4.__version__,
        platform.python_version(),
    )

    ds = insert_history(ds, histtext)

    return ds


def set_var_dtype(ds, var, dtype="float32"):
    ds[var].encoding["dtype"] = dtype

    return ds


def open_time_2d_dataset(filename):
    # need to drop 'time' variable because of xarray limitations related
    # to coordinates and variables with the same name, otherwise it raises a
    # MissingDimensionsError
    # Check if CF or not, and return the correct dataset
    with xr.open_dataset(filename, decode_times=False, drop_variables="time") as ds:
        if is_cf(ds):
            iscf = True
        else:
            iscf = False

    if iscf:
        return xr.open_dataset(filename)
    else:
        return xr.open_dataset(filename, decode_times=False, drop_variables="time")


def epic_to_cf_time(ds):
    if "time_cf" in ds:
        ds["time"] = ds["time_cf"]
    else:
        ds["time"] = epic_to_datetime(ds["time"].values, ds["time2"].values)

    for v in ["time_cf", "time2", "epic_time", "epic_time2"]:
        if v in ds:
            ds = ds.drop(v)
    return xr.decode_cf(ds, decode_times=True)


def epic_to_datetime(time, time2):
    thedate = pd.to_datetime(time - 0.5, origin="julian", unit="D")
    thetime = pd.to_timedelta(time2, unit="ms")
    return thedate + thetime


def make_jd(time):
    return time.to_julian_date().values + 0.5


def make_epic_time(jd):
    epic_time = np.floor(jd)
    # make sure they are all integers, and then cast as such
    if np.all(np.mod(epic_time, 1) == 0):
        epic_time = epic_time.astype(np.int32)
    else:
        warnings.warn(
            "Not all EPIC time values are integers; "
            "this will cause problems with time and time2"
        )

    return epic_time


def make_epic_time2(jd):
    # TODO: Hopefully this is correct... roundoff errors on big numbers...
    return np.round((jd - np.floor(jd)) * 86400000).astype(np.int32)


def create_epic_times(ds, waves=False):
    jd = make_jd(ds["time"].to_dataframe().index)

    ds["epic_time"] = xr.DataArray(make_epic_time(jd), dims="time")
    ds["epic_time"].encoding["_FillValue"] = None

    ds["epic_time2"] = xr.DataArray(make_epic_time2(jd), dims="time")
    ds["epic_time2"].encoding["_FillValue"] = None

    return ds


def check_update_attrs(ds, key, value):
    """Update attr and raise warning if attr already exists and is different from replacement value"""
    if key in ds.attrs and ds.attrs[key] != value:
        warnings.warn(f"attrs collision. Replacing {ds.attrs[key]=} with '{value}'.")

    ds.attrs[key] = value

    return ds


def add_start_stop_time(ds):
    """Add start_time and stop_time attrs"""

    ds = check_update_attrs(ds, "start_time", ds["time"][0].values.astype(str))
    ds = check_update_attrs(ds, "stop_time", ds["time"][-1].values.astype(str))

    return ds


# def add_lat_lon(ds, var):
#     """Add lat and lon dimensions"""
#
#     ds[var] = xr.concat([ds[var]], dim=ds["lon"])
#     ds[var] = xr.concat([ds[var]], dim=ds["lat"])
#
#     # Reorder so lat, lon are at the end.
#     dims = [d for d in ds[var].dims if (d != "lon") and (d != "lat")]
#     dims.extend(["lat", "lon"])
#     dims = tuple(dims)
#
#     ds[var] = ds[var].transpose(*dims)
#
#     return ds


def ds_add_lat_lon(ds):
    ds["latitude"] = xr.DataArray(
        [ds.attrs["latitude"]],
        dims=("latitude"),
        name="latitude",
        attrs={
            "units": "degree_north",
            "axis": "Y",
            "long_name": "Latitude",
            "standard_name": "latitude",
            "epic_code": 500,
        },
    )

    ds["longitude"] = xr.DataArray(
        [ds.attrs["longitude"]],
        dims=("longitude"),
        name="longitude",
        attrs={
            "units": "degree_east",
            "axis": "X",
            "long_name": "Longitude",
            "standard_name": "longitude",
            "epic_code": 502,
        },
    )

    return ds


def shift_time(ds, timeshift, apply_clock_error=True, apply_clock_drift=True):
    """Shift time to middle of burst and apply clock error/offset and drift correction"""

    if timeshift != 0:
        # back up attrs as these are lost in the process
        attrsbak = ds["time"].attrs
        # shift times to center of ensemble
        ds["time"] = ds["time"] + np.timedelta64(int(timeshift), "s")
        ds["time"].attrs = attrsbak

        if not timeshift.is_integer():
            warnings.warn(
                "time offset of %.3f s was adjusted to %.f s for shifting time"
                % (timeshift, int(timeshift))
            )

        histtext = "Time shifted to middle of burst by {} s.".format(
            int(timeshift),
        )

        insert_history(ds, histtext)

    if "ClockError" in ds.attrs and apply_clock_error:
        if ds.attrs["ClockError"] != 0:
            # back up attrs as these are lost in the process
            attrsbak = ds["time"].attrs
            # note negative on ds.attrs['ClockError']
            ds["time"] = ds["time"] + np.timedelta64(-ds.attrs["ClockError"], "s")
            ds["time"].attrs = attrsbak

            histtext = "Time shifted by {:d} s from ClockError.".format(
                -ds.attrs["ClockError"],
            )

            insert_history(ds, histtext)

    if "ClockDrift" in ds.attrs and apply_clock_drift:
        if ds.attrs["ClockDrift"] != 0:
            # back up attrs as these are lost in the process
            attrsbak = ds["time"].attrs
            # note negative on ds.attrs['ClockDrift']
            ds["time"] = ds["time"] + pd.TimedeltaIndex(
                np.linspace(0, -ds.attrs["ClockDrift"], len(ds["time"])), "s"
            )

            ds["time"] = ds.time.dt.round("1s")
            ds["time"].attrs = attrsbak

            histtext = "Time linearly interpolated by {} s using ClockDrift and rounded to the nearest second.\n".format(
                -ds.attrs["ClockDrift"],
            )

            insert_history(ds, histtext)

    return ds


def create_water_depth_var(ds, psh="initial_instrument_height"):
    press = None

    if "Pressure_ac" in ds:
        press = "Pressure_ac"
    elif "P_1ac" in ds:
        press = "P_1ac"
    elif "Pressure" in ds:
        press = "Pressure"
    elif "P_1" in ds:
        press = "P_1"

    if press is not None:
        if "sample" in ds.dims:
            ds["water_depth"] = xr.DataArray(
                ds[press].squeeze().mean(dim="sample") + ds.attrs[psh]
            )

        else:
            ds["water_depth"] = xr.DataArray(ds[press].squeeze() + ds.attrs[psh])

        ds["water_depth"].attrs["long_name"] = "Total water depth"
        ds["water_depth"].attrs["units"] = "m"
        ds["water_depth"].attrs["standard_name"] = "sea_floor_depth_below_sea_surface"
        ds["water_depth"].attrs["epic_code"] = 3

        histtext = f"Create water_depth variable using {press} and {psh} attribute"
        ds = insert_history(ds, histtext)

    return ds


# def create_water_depth(ds):
#     """Create WATER_DEPTH attribute"""
#
#     press = None
#
#     if "Pressure_ac" in ds:
#         press = "Pressure_ac"
#     elif "P_1ac" in ds:
#         press = "P_1ac"
#     elif "Pressure" in ds:
#         press = "Pressure"
#     elif "P_1" in ds:
#         press = "P_1"
#
#     if "sample" in ds.dims:
#         dims = ("time", "sample")
#     else:
#         dims = "time"
#
#     if "initial_instrument_height" in ds.attrs:
#         if press:
#             ds.attrs["nominal_instrument_depth"] = (
#                 ds[press].squeeze().mean(dim=dims).values
#             )
#             # ds['water_depth'] = ds.attrs['nominal_instrument_depth']
#             wdepth = (
#                 ds.attrs["nominal_instrument_depth"]
#                 + ds.attrs["initial_instrument_height"]
#             )
#             if "ac" in press:
#                 ds.attrs["WATER_DEPTH_source"] = (
#                     "water depth = MSL from "
#                     "pressure sensor, "
#                     "atmospherically corrected"
#                 )
#             else:
#                 ds.attrs["WATER_DEPTH_source"] = (
#                     "water depth = MSL from " "pressure sensor"
#                 )
#             ds.attrs["WATER_DEPTH_datum"] = "MSL"
#         else:
#             wdepth = ds.attrs["WATER_DEPTH"]
#             ds.attrs["nominal_instrument_depth"] = (
#                 ds.attrs["WATER_DEPTH"] - ds.attrs["initial_instrument_height"]
#             )
#         # ds['Depth'] = ds.attrs['nominal_instrument_depth']
#         # TODO: why is this being redefined here? Seems redundant
#         ds.attrs["WATER_DEPTH"] = wdepth
#
#     elif "nominal_instrument_depth" in ds.attrs:
#         ds.attrs["initial_instrument_height"] = (
#             ds.attrs["WATER_DEPTH"] - ds.attrs["nominal_instrument_depth"]
#         )
#         # ds['water_depth'] = ds.attrs['nominal_instrument_depth']
#
#     if "initial_instrument_height" not in ds.attrs:
#         # TODO: do we really want to set to zero?
#         ds.attrs["initial_instrument_height"] = 0
#
#     return ds


def create_nominal_instrument_depth(ds):
    if "nominal_instrument_depth" not in ds.attrs:
        ds.attrs["nominal_instrument_depth"] = (
            ds.attrs["WATER_DEPTH"] - ds.attrs["initial_instrument_height"]
        )

    return ds


def create_z(ds):
    # create z and depth coordinate variables
    iih = ds.attrs["initial_instrument_height"]
    if "bindist" in ds:
        bd = ds["bindist"].values

    if "pressure_sensor_height" in ds.attrs:
        psh = ds.attrs["pressure_sensor_height"]
    else:
        psh = None

    if "NAVD88_ref" in ds.attrs:
        hagd = ds.attrs["NAVD88_ref"]
        gdn = "NAVD88"
    elif "height_above_geopotential_datum" in ds.attrs:
        hagd = ds.attrs["height_above_geopotential_datum"]
        gdn = ds.attrs["geopotential_datum_name"]
    else:
        hagd = 0
        gdn = "sea bed"

    if "bindist" in ds:
        if ds.attrs["orientation"].upper() == "DOWN":
            ds["z"] = xr.DataArray(hagd + iih - bd, dims="z")
        elif ds.attrs["orientation"].upper() == "UP":
            ds["z"] = xr.DataArray(hagd + iih + bd, dims="z")
    else:
        ds["z"] = xr.DataArray([hagd + iih], dims="z")

    if gdn != "sea bed":
        ds["z"].attrs["geopotential_datum_name"] = gdn
    ds["z"].attrs["long_name"] = f"height relative to {gdn}"

    ds["z"].attrs["positive"] = "up"
    ds["z"].attrs["axis"] = "Z"
    ds["z"].attrs["units"] = "m"
    ds["z"].attrs["standard_name"] = "height"

    if psh is not None:
        # Create zsen for sensor data
        ds["zsen"] = xr.DataArray([hagd + psh], dims="zsen")
        ds["zsen"].attrs["positive"] = "up"
        ds["zsen"].attrs["units"] = "m"
        ds["zsen"].attrs["standard_name"] = "height"
        if gdn != "sea bed":
            ds["zsen"].attrs["geopotential_datum_name"] = gdn
        ds["zsen"].attrs["long_name"] = f"pressure sensor height relative to {gdn}"

    # find depth dimension values
    # check for sea floor depth standard names
    # this creates a depth variable for the instrument measurement location
    # if there are multiple attributes in the config file, it will prefer ones higher in the list below
    sfds = [
        "sea_floor_depth_below_geoid",
        "sea_floor_depth_below_geopotential_datum",
        "sea_floor_depth_below_mean_sea_level",
        "sea_floor_depth_below_reference_ellipsoid",
        "sea_floor_depth_below_sea_surface",
        "WATER_DEPTH",
    ]
    depvar = None
    for name in sfds:
        if name in ds.attrs:
            if "sea_floor" in name:
                longname = name.replace("sea_floor_", "").replace("_", " ")
            else:
                longname = f"depth below {name} var"
            if "bindist" in ds:
                if ds.attrs["orientation"].upper() == "DOWN":
                    depvar = ds.attrs[name] - iih + bd
                elif ds.attrs["orientation"].upper() == "UP":
                    depvar = ds.attrs[name] - iih - bd
            else:
                depvar = [ds.attrs[name] - iih]

            if psh is not None:
                depsenvar = [ds.attrs[name] - psh]

        if depvar is not None:
            break

    if depvar is None:
        # only if we don't already have a sea_floor_depth_below_<X> specified above
        if "D_3" in ds or "P_1ac" in ds:
            # here we es
            if "D_3" in ds:
                v = "D_3"
            elif "P_1ac" in ds:
                v = "P_1ac"

            if "bindist" in ds:
                if ds.attrs["orientation"].upper() == "DOWN":
                    depvar = np.nanmean(ds[v]) + bd
                elif ds.attrs["orientation"].upper() == "UP":
                    depvar = np.nanmean(ds[v]) - bd
            else:
                depvar = [np.nanmean(ds[v])]

            # name = "sea_floor_depth_below_mean_sea_level"
            longname = "depth below mean sea level from data"

            if psh is not None:
                depsenvar = [np.nanmean(ds[v])]

    ds["depth"] = xr.DataArray(depvar, dims="depth")

    ds["depth"].attrs["positive"] = "down"
    ds["depth"].attrs["units"] = "m"
    ds["depth"].attrs["standard_name"] = "depth"
    ds["depth"].attrs["long_name"] = longname

    if psh is not None:
        # Create depthsen for sensor data
        ds["depthsen"] = xr.DataArray(depsenvar, dims="depthsen")
        ds["depthsen"].attrs["positive"] = "down"
        ds["depthsen"].attrs["units"] = "m"
        ds["depthsen"].attrs["standard_name"] = "depth"
        ds["depthsen"].attrs["long_name"] = f"pressure sensor {longname}"

    return ds


# def add_z_if_no_pressure(ds, var):
#     # no_p = no pressure sensor. also use for exo
#     attrsbak = ds["z"].attrs
#     ds[var] = ds[var].expand_dims("z")
#     # reorder so z at end
#     dims = [d for d in ds[var].dims if (d != "z")]
#     dims.extend(["z"])
#     dims = tuple(dims)
#     ds[var] = ds[var].transpose(*dims)
#     ds["z"].attrs = attrsbak
#
#     return ds


# def no_p_create_depth(ds):
#     # no_p = no pressure sensor. also use for exo
#     if "NAVD88_ref" in ds.attrs:
#         ds["depth"] = xr.DataArray(
#             [-ds.attrs["NAVD88_ref"] - ds.attrs["initial_instrument_height"]],
#             dims="depth",
#         )
#         ds["depth"].attrs["VERT_DATUM"] = "NAVD88"
#         ds["depth"].attrs["NOTE"] = (
#             "Computed as platform depth "
#             "[m NAVD88] minus "
#             "initial_instrument_height"
#         )
#     else:
#         ds["depth"] = xr.DataArray(
#             [ds.attrs["WATER_DEPTH"] - ds.attrs["initial_instrument_height"]],
#             dims="depth",
#         )
#         ds["depth"].attrs["NOTE"] = (
#             "Computed as WATER_DEPTH minus " "initial_instrument_height"
#         )
#
#     ds["depth"].attrs["positive"] = "down"
#     ds["depth"].attrs["axis"] = "Z"
#     ds["depth"].attrs["units"] = "m"
#     ds["depth"].attrs["epic_code"] = 3
#     ds["depth"].encoding["_FillValue"] = None
#
#     return ds


# def no_p_add_depth(ds, var):
#     # no_p = no pressure sensor. also use for exo
#     ds[var] = xr.concat([ds[var]], dim=ds["depth"])
#
#     # Reorder so lat, lon are at the end.
#     dims = [d for d in ds[var].dims if (d != "depth")]
#     dims.extend(["depth"])
#     dims = tuple(dims)
#
#     ds[var] = ds[var].transpose(*dims)
#
#     return ds


def insert_note(ds, var, notetxt):
    toinsert = "{}: {}\n".format(
        datetime.datetime.now(datetime.timezone.utc).isoformat(), notetxt
    )

    print(f"{var}: {toinsert.rstrip()}")

    if "note" in ds[var].attrs:
        ds[var].attrs["note"] = ds[var].attrs["note"] + toinsert
    else:
        ds[var].attrs["note"] = toinsert

    return ds


def add_delta_t(ds):
    deltat = ((ds["time"][1] - ds["time"][0]) / np.timedelta64(1, "s")).item()
    if not deltat.is_integer():
        warnings.warn("DELTA_T is not an integer; casting as int in attrs")

    ds.attrs["DELTA_T"] = int(round(deltat))

    return ds


def atmos_correct(ds, atmpres):
    met = xr.load_dataset(atmpres)
    # need to save attrs before the subtraction, otherwise they are lost
    attrs = ds["P_1"].attrs
    ds["P_1ac"] = ds["P_1"] - met["atmpres"] - met["atmpres"].offset
    print(
        f"Atmospherically correcting using time-series from {atmpres} and offset of {met['atmpres'].offset}"
    )
    ds["P_1ac"].attrs = attrs

    ds.attrs["atmospheric_pressure_correction_file"] = atmpres
    ds.attrs["atmospheric_pressure_correction_offset_applied"] = met["atmpres"].attrs[
        "offset"
    ]
    if "comment" in met["atmpres"].attrs:
        ds.attrs["atmospheric_pressure_correction_comment"] = met["atmpres"].attrs[
            "comment"
        ]

    # Is reindexing with a tolerance still necessary? Keeping old code around to check
    # ds["P_1ac"] = (
    #     ds["Pressure"]
    #     - met["atmpres"].reindex_like(
    #         ds["Pressure"], method="nearest", tolerance="10min"
    #     )
    #     - met["atmpres"].offset
    # )

    return ds


def read_samplingrates_burst(ds, conn):
    """
    Reads in sample information from RBR instrument in burst mode
    """
    # Get samples per burst
    try:
        # this seems to be used on older-style databases;
        # throws error on newer files
        samplingcount = conn.execute("select samplingcount from schedules").fetchall()[
            0
        ][0]
    except sqlite3.OperationalError:
        samplingcount = conn.execute("select samplingcount from wave").fetchall()[0][0]
    ds.attrs["samples_per_burst"] = samplingcount

    # Get sampling interval
    try:
        samplingperiod = conn.execute(
            "select samplingperiod from schedules"
        ).fetchall()[0][0]
    except sqlite3.OperationalError:
        samplingperiod = conn.execute("select samplingperiod from wave").fetchall()[0][
            0
        ]
    ds.attrs["sample_interval"] = samplingperiod / 1000

    # Get Repetition interval of sampling
    try:
        repetitionperiod = conn.execute(
            "select repetitionperiod from schedules"
        ).fetchall()[0][0]
    except sqlite3.OperationalError:
        repetitionperiod = conn.execute("select repetitionperiod from wave").fetchall()[
            0
        ][0]

    # Convert to seconds
    ds.attrs["burst_interval"] = repetitionperiod / 1000

    # Length of bursts in data points
    ds.attrs["burst_length"] = ds.attrs["samples_per_burst"]

    return ds


def read_samplingrates_continuous(ds, conn):
    """
    Reads in sample information from RBR instrument in continuous mode
    """
    try:
        samplingperiod = conn.execute(
            "select samplingperiod from schedules"
        ).fetchall()[0][0]
    except sqlite3.OperationalError:
        samplingperiod = conn.execute(
            "select samplingperiod from continuous"
        ).fetchall()[0][0]

    samplingperiod = samplingperiod / 1000  # convert from ms to sec
    samplingrate = 1 / samplingperiod  # convert to rate [Hz]

    # Set sampling period, [sec]
    ds.attrs["sample_interval"] = samplingperiod

    # Set samples per burst
    samplingcount = ds.attrs["wave_interval"] * samplingrate
    ds.attrs["samples_per_burst"] = round(samplingcount)

    # Set burst interval, [sec], USER DEFINED in instrument attr
    ds.attrs["burst_interval"] = ds.attrs["wave_interval"]

    # Set sample interval
    ds.attrs["burst_length"] = ds.attrs["samples_per_burst"]

    return ds


def salinity_from_spcon(spcon):
    """
    Derive salinity from specific conductance following
    Simplified conversions between specific conductance and salinity units for use with data from monitoring stations
    Interagency Ecological Program Newsletter
    By: Laurence E. Schemel
    https://pubs.er.usgs.gov/publication/70174311
    spcon must be in uS/cm
    """

    R = spcon / 53087
    K1 = 0.0120
    K2 = -0.2174
    K3 = 25.3283
    K4 = 13.7714
    K5 = -6.4788
    K6 = 2.5842

    return K1 + K2 * R**0.5 + K3 * R + K4 * R ** (3 / 2) + K5 * R**2 + K6 * R ** (5 / 2)


def spcon_from_salinity(S):
    """
    Derive specific conductance from salinity following
    Simplified conversions between specific conductance and salinity units for use with data from monitoring stations
    Interagency Ecological Program Newsletter
    By: Laurence E. Schemel
    https://pubs.er.usgs.gov/publication/70174311
    Returns spcon in uS/cm
    """

    J1 = -16.072
    J2 = 4.1495
    J3 = -0.5345
    J4 = 0.0261

    return S / 35 * 53087 + S * (S - 35) * (
        J1 + (J2 * S ** (1 / 2)) + (J3 * S) + (J4 * S ** (3 / 2))
    )


def check_fits_in_int32(ds, var):
    if np.nanmax(np.abs(ds[var])) > (2**31 - 1):
        warnings.warn(
            f"32-bit integer overflow on {var}; setting encoding to i4 will fail"
        )
        return False
    else:
        return True


def check_time_fits_in_int32(ds, var):
    num, units, calendar = xr.coding.times.encode_cf_datetime(ds[var])
    if np.nanmax(num) > (2**31 - 1):
        warnings.warn(
            f"32-bit integer overflow on {var}; setting encoding to i4 will fail"
        )
        return False
    else:
        return True


def check_time_encoding(ds):
    """ensure we don't set dtypes uint for CF compliance"""

    if "units" in ds["time"].encoding:
        ds["time"].encoding.pop("units")

    if check_time_fits_in_int32(ds, "time"):
        ds["time"].encoding["dtype"] = "i4"

    else:
        print("time variable will not fit in int32; casting to double")
        ds["time"].encoding["dtype"] = "double"

    return ds


def check_valid_globalatts_metadata(metadata):
    for k in ["WATER_DEPTH", "latitude", "longitude", "MOORING"]:
        if k not in metadata:
            raise KeyError(
                f"{k} must be defined, most likely in global attributes file"
            )


def read_globalatts(fname):
    """
    Read global attributes file (glob_attxxxx.txt) and create metadata
    structure
    """

    metadata = {}

    with open(fname) as csvfile:
        a = csv.reader(csvfile, delimiter=";")

        for row in a:
            if row[0].strip() == "MOORING":
                metadata[row[0].strip()] = str(row[1].strip())
            else:
                metadata[row[0].strip()] = str2num(row[1].strip())

        return metadata


def str2num(s):
    """
    Convert string to float if possible
    """

    try:
        float(s)
        return float(s)
    except ValueError:
        return s


def loadmat(filename):
    """
    this function should be called instead of direct spio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects

    from: `StackOverflow <https://stackoverflow.com/q/7008608>`_
    """
    data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)


def _check_keys(dic):
    """
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries
    """
    for key in dic:
        if isinstance(dic[key], spio.matlab.mio5_params.mat_struct):
            dic[key] = _todict(dic[key])
    return dic


def _todict(matobj):
    """
    A recursive function which constructs from matobjects nested dictionaries
    """
    dic = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, spio.matlab.mio5_params.mat_struct):
            dic[strg] = _todict(elem)
        else:
            dic[strg] = elem
    return dic


def create_water_level_var(ds, var="P_1ac", vdim="z"):
    """
    Create water level variable from NAVD88 sensor height
    """

    if (
        var in ds.data_vars
        and "geopotential_datum_name" in ds[vdim].attrs
        and ds[vdim].attrs["geopotential_datum_name"] == "NAVD88"
    ):
        if "sample" in ds.dims:
            ds["water_level"] = xr.DataArray(
                ds[var].squeeze().mean(dim="sample") + ds[vdim].values
            )
        else:
            ds["water_level"] = ds[var] + ds[vdim].values

        ds["water_level"].attrs["long_name"] = "Water level NAVD88"
        ds["water_level"].attrs["units"] = "m"
        ds["water_level"].attrs[
            "standard_name"
        ] = "sea_surface_height_above_geopotential_datum"
        ds["water_level"].attrs["geopotential_datum_name"] = "NAVD88"
        ds["water_level"].attrs["note"] = f"Water level calculated using {var} + {vdim}"

        histtext = f"Create water_level variable relative to NAVD88 using {var} and {vdim} vertical dimension"
        ds = insert_history(ds, histtext)

    else:
        print(
            f"Cannot create water_level variable without {var} and height_above_geopotential_datum relative to NAVD88 in global attributes file."
        )
    return ds


def create_filtered_water_level_var(ds):
    """
    Create 4th order lowpass butterworth filtered water level with 6 min cutoff
    """

    if "filtered_wl" in ds.attrs and ds.attrs["filtered_wl"].lower() == "true":

        var = "water_level"
        cutfreq = 1 / 360  # 6 min cutoff
        ftype = "lowpass"
        ford = 4

        if "sample_rate" in ds.attrs:
            sr = ds.attrs["sample_rate"]
        elif "sample_interval" in ds.attrs:
            sr = 1 / ds.attrs["sample_interval"]
        else:
            raise ValueError(
                "Cannot create filtered_water_level without sample_rate or sample _interval in global attributes"
            )

        filtered_wl = filter.butter_filt(ds[var].values, sr, cutfreq, ftype, ford)

        ds["water_level_filt"] = xr.DataArray(filtered_wl, dims="time")

        ds["water_level_filt"].attrs["long_name"] = "Filtered water level NAVD88"
        ds["water_level_filt"].attrs["units"] = "m"
        ds["water_level_filt"].attrs[
            "standard_name"
        ] = "sea_surface_height_above_geopotential_datum"
        ds["water_level_filt"].attrs["geopotential_datum_name"] = "NAVD88"
        ds["water_level_filt"].attrs[
            "note"
        ] = "4th order lowpass butterworth filter with 6 min cutoff"

    return ds


def rename_diwasp_wave_vars(ds):
    """
    Rename DIWASP variables to standard EPIC names used in stglib wave statistics
    """
    # rename wave vars to user specified convention
    if ds.attrs["diwasp_names"].lower() == "epic":
        varnames = {
            "diwasp_frequency": "frequency",
            "diwasp_direction": "direction",
            "diwasp_hs": "wh_4061",
            "diwasp_tp": "wp_peak",
            "diwasp_tm": "wp_4060",
            "diwasp_dtp": "wvdir",
            "diwasp_dp": "dwvdir",
            "diwasp_dm": "wd_4062",
            "diwasp_pspec": "pspec",
            "diwasp_astspec": "sspec",
            "diwasp_vspec": "vspec",
            "diwasp_dspec": "dspec",
        }

        # check to make sure they exist before trying to rename
        newvars = {}
        for k in varnames:
            if k in ds:
                newvars[k] = varnames[k]
    else:
        return ds

    return ds.rename(newvars)


def rename_diwasp_fspec(diwasp):
    """Rename diwasp_fspec using first input datatype"""

    if "diwasp_inputs" in diwasp.attrs:
        inputs = diwasp.attrs["diwasp_inputs"]
        newname = {}
        if inputs[0] == "elev":
            newname = {"diwasp_fspec": "diwasp_astspec"}
        elif inputs[0] == "pres":
            newname = {"diwasp_fspec": "diwasp_pspec"}
        elif inputs[0] in ["velx", "vely", "velz", "radial"]:
            newname = {"diwasp_fspec": "diwasp_vspec"}
        else:
            return diwasp
    else:
        return diwasp

    return diwasp.rename(newname)


def clip_ds_prf(ds):
    """
    Clip profile data set by user specified bindist range - looks for vertical coordinates ["z", "depth", "bindist", "bins" ] to clip
    """

    # Find valid coordinate dimensions for dataset
    coords = []
    for k in ds.coords:
        coords.append(k)

    if "good_bindist" in ds.attrs:
        if "bindist" in coords:
            idx1 = int(
                ds["bindist"]
                .where(ds.bindist > ds.attrs["good_bindist"][0])
                .argmin()
                .values
            )
            idx2 = int(
                ds["bindist"]
                .where(ds.bindist > ds.attrs["good_bindist"][1])
                .argmin()
                .values
            )

            for k in ["z", "bindist", "depth", "bins"]:
                if k in ds.coords:
                    ds = ds.isel({k: slice(idx1, idx2)})

            histtext = f"Profile data clipped by bin distance range {ds.attrs['good_bindist'][0]} m to {ds.attrs['good_bindist'][1]} m"
            ds = insert_history(ds, histtext)

        else:

            print(
                f"Profile data not clipped using 'good_bindist', 'bindist' is not a valid coordinate dimension {coords}"
            )

    else:

        print("Profile data not clipped; no values specified in metadata")

    return ds


def cart2pol(x, y):
    """
    Simple cartesian to polar coordinate transformation
    returns:
    rho - radius in same units as x and y
    phi - angle in radians
    """
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return (rho, phi)


def pol2cart(rho, phi):
    """
    Simple polar to cartesian coordinate transformation
    inputs:
    rho - radius in same units as returned values for x and y
    phi - angle in radians
    """
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return (x, y)


def make_vector_average_vars(ds, data_vars=None, dim=None):
    """
    Find vector averages for specified variables, along specified dimension
    """
    dsVA = xr.Dataset()
    if data_vars is not None and dim is not None:
        for var in data_vars:
            attrsbak = []
            if var in ds.data_vars:
                attrsbak = ds[var].attrs
                R = ds[var]
                Rx, Ry = pol2cart(1, np.radians(R))
                rho, phi = cart2pol(Rx.mean(dim=dim), Ry.mean(dim=dim))
                dsVA[var] = np.degrees(phi)
                dsVA[var].attrs = attrsbak

    return dsVA
