import gsw
import numpy as np
import pandas as pd
import xarray as xr

from .aqd import aqdutils
from .core import qaqc, utils


def log_to_cdf(metadata):
    basefile = metadata["basefile"]

    # if ea400
    if metadata["instrument_type"] == "ea":
        instmeta = read_ea_instmet(basefile + ".log")
        ds = load_ea_point(basefile + ".log", instmeta)
        ds = load_ea_profile(ds, basefile + ".log", instmeta)

    # if aa400
    elif metadata["instrument_type"] == "aa":
        instmeta = read_aa_instmet(basefile + ".log")
        ds = load_aa_point(
            basefile + ".log",
            instmeta,
            skiprows=metadata["skiprows"],
            skipfooter=metadata["skipfooter"],
        )
        metadata.pop("skiprows")
        metadata.pop("skipfooter")

    else:
        raise ValueError(
            "instrument_type in config file, {:s}, is invalid".format(
                metadata["instrument_type"]
            )
        )

    prefix = metadata["instrument_type"].upper()

    metadata["instmeta"] = instmeta

    ds = utils.write_metadata(ds, metadata)
    ds = utils.ensure_cf(ds)
    ds = utils.shift_time(ds, 0)

    if ds.attrs[f"{prefix}Pulses_in_series_num"] > 1:
        ds.attrs.update({"sample_mode": "BURST"})
    else:
        ds.attrs.update({"sample_mode": "CONTINUOUS"})

    # configure file
    print("Configuring .cdf file")

    cdf_filename = ds.attrs["filename"] + "-raw.cdf"

    ds.to_netcdf(
        cdf_filename,
        unlimited_dims=["time"],
    )

    print("Finished writing data to %s" % cdf_filename)

    return ds


def cdf_to_nc(cdf_filename):
    """
    Load a raw .cdf file and generate a processed .nc file
    """

    # Load raw .cdf data
    ds = xr.load_dataset(cdf_filename)

    if "bins" in ds:
        ds = calc_bin_height(ds)  # calculate bin height
        # calculate corrected bin height (on NAVD88 datum) with adjusted sound speed
        ds = calc_cor_bin_height(ds)
        # Trim data when instrument is out of water during a deployment and extra bins if good_bins specified
        ds = trim_alt(ds)
    else:
        ds = trim_alt(ds, data_vars=["Altitude_m", "AmplitudeFS", "Temperature_C"])

    # calculate corrected altitude (distance to bed/b_range) with adjusted sound speed
    ds = calc_cor_brange(ds)

    # calculate seabed elevation referenced to datum
    ds = calc_boundary_elev(ds)

    ds = utils.create_z(ds)

    # swap bin dim with bin_height
    ds = aqdutils.ds_swap_dims(ds)  # swap vert dim to z or user specified in vert_dim

    # rename variables
    ds = ds_rename_vars(ds)

    # drop some vars after renaming
    for k in [
        "bins",
        "Ping",
        "Ping_num_in_series",
        "Altitude_m",
        "Battery_mV",
    ]:
        if k in ds:
            ds = ds.drop_vars(k)

    # add attributes to each variable
    ds = ds_add_attrs(ds)

    # Call utils
    ds = utils.clip_ds(ds)
    ds = utils.add_min_max(ds)
    ds = utils.ds_add_lat_lon(ds)
    ds = utils.ds_coord_no_fillvalue(ds)
    ds = utils.add_start_stop_time(ds)
    ds = utils.check_time_encoding(ds)
    ds = utils.add_delta_t(ds)

    # Write to .nc file
    print("Writing cleaned/trimmed burst data and averaged burst data to .nc file")
    nc_filename = ds.attrs["filename"] + "b.nc"

    ds.to_netcdf(nc_filename, unlimited_dims=["time"])

    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])

    print("Done writing netCDF file", nc_filename)

    # Average burst and write to -a.nc file
    ds = average_burst(ds)

    # QAQC
    ds = qaqc.call_qaqc(ds)

    # Call utils
    ds = utils.add_min_max(ds)

    nc_filename = ds.attrs["filename"] + "b-a.nc"

    ds.to_netcdf(
        nc_filename,
        unlimited_dims=["time"],
    )

    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])

    print("Done writing burst averaged netCDF file", nc_filename)

    return ds


def read_ea_instmet(basefile):
    with open(basefile) as f:
        instmeta = {}
        row = ""
        while "##DataStart" not in row:
            row = f.readline().rstrip()
            dat = row.split()
            if "#DeviceID" in row:
                instmeta["EADeviceID"] = row[10:]
                instmeta["serial_number"] = dat[1]
            elif "#NSamples" in row:
                instmeta["EABin_count"] = int(dat[1])
            elif "#Resolution,m" in row:
                instmeta["EABin_size_m"] = float(dat[1])
            elif "#SoundSpeed,mps" in row:
                instmeta["EASoundSpeed_mps"] = float(dat[1])
            elif "#Tx_Frequency,Hz" in row:
                instmeta["EATx_Frequency_Hz"] = float(dat[1])
            elif "#Range,m" in row:
                instmeta["EARange_m"] = float(dat[1])
            elif "#Pulse period,sec" in row:
                instmeta["EAPulse_period_sec"] = float(dat[2])
            elif "#Pulses in series,num" in row:
                instmeta["EAPulses_in_series_num"] = int(dat[3])
            elif "#Interval between series,sec" in row:
                instmeta["EAInterval_between_series_sec"] = float(dat[3])
            elif "#Threshold,%" in row:
                instmeta["EAThreshold_percent"] = int(dat[1])
            elif "#Offset,m" in row:
                instmeta["EAOffset_m"] = float(dat[1])
            elif "#Deadzone,m" in row:
                instmeta["EADeadzone_m"] = float(dat[1])
            elif "#PulseLength,uks" in row:
                instmeta["EAPulseLength_microsec"] = float(dat[1])
            elif "#TVG_Gain,dB" in row:
                instmeta["EATVG_Gain_dB"] = float(dat[1])
            elif "#TVG_Slope,dB/km" in row:
                instmeta["EATVG_Slope_dBkm"] = float(dat[1])
            elif "#TVG_Mode" in row:
                instmeta["EATVG_Mode"] = int(dat[1])
            elif "#OutputMode" in row:
                instmeta["EAOutputMode"] = int(dat[1])

    return instmeta


def load_ea_point(basefile, instmeta):
    with open(basefile) as f:
        data = f.read().splitlines()

        point = {}
        point["TimeLocal"] = []
        point["TimeUTC"] = []
        point["Ping"] = []
        point["Ping_num_in_series"] = []
        point["Altitude_m"] = []
        point["Temperature_C"] = []
        point["Pitch_deg"] = []
        point["Roll_deg"] = []

        for row in data:
            dat = row.split()
            if "#TimeLocal" in row:
                point["TimeLocal"].append(dat[1] + " " + dat[2])
            elif "#TimeUTC" in row:
                point["TimeUTC"].append(dat[1] + " " + dat[2])
            elif "#Ping  " in row:
                point["Ping"].append(float(dat[1]))
            elif "#Ping num in series" in row:
                point["Ping_num_in_series"].append(float(dat[4]))
            elif "#Altitude,m" in row:
                point["Altitude_m"].append(float(dat[1]))
            elif "#Temperature" in row:
                point["Temperature_C"].append(float(dat[1]))
            elif "#Pitch,deg" in row:
                point["Pitch_deg"].append(float(dat[1]))
            elif "#Roll,deg" in row:
                point["Roll_deg"].append(float(dat[1]))

    # point["TimeUTC"] = pd.to_datetime(point["TimeUTC"])
    point["TimeUTC"] = pd.to_datetime(point["TimeUTC"]).to_numpy()

    samples = instmeta["EAPulses_in_series_num"]
    n = instmeta["EABin_count"]

    for k in point:
        point[k] = np.array(point[k])
        point[k] = point[k].reshape((-1, samples))

    time = point["TimeUTC"][:, 0]

    ds = xr.Dataset()
    ds["time"] = xr.DataArray(time, dims="time")
    ds["sample"] = xr.DataArray(np.arange(0, samples), dims="sample")
    ds["bins"] = xr.DataArray(np.arange(n), dims="bins")

    # add variables to xarray
    for k in point:
        if "Time" not in k:
            ds[k] = xr.DataArray(point[k], dims=("time", "sample"))

    return ds


def load_ea_profile(ds, basefile, instmeta):
    profile = []
    with open(basefile) as f:  # read in profile echo data
        for row in f:
            if row.rstrip() == "##DataStart":
                for row in f:
                    if row.rstrip() == "##DataEnd":
                        break
                    profile.append(float(row))

    # reshape profile data
    samples = instmeta["EAPulses_in_series_num"]
    n = instmeta["EABin_count"]
    profile = np.array(profile)
    profile = profile.reshape((-1, samples, n))
    ds["Counts"] = xr.DataArray(
        profile,
        dims=("time", "sample", "bins"),
        coords=[ds["time"], ds["sample"], ds["bins"]],
    )

    return ds


def ds_rename_vars(ds):
    varnames = {
        "Temperature_C": "Tx_1211",
        "Pitch_deg": "Ptch_1216",
        "Roll_deg": "Roll_1217",
        "Counts": "AGC_1202",
        "AmplitudeFS": "AMP_723",
    }

    # check to make sure they exist before trying to rename
    newvars = {}
    for k in varnames:
        if k in ds:
            newvars[k] = varnames[k]

    return ds.rename(newvars)


def ds_add_attrs(ds):

    ds["time"].attrs.update(
        {"standard_name": "time", "axis": "T", "long_name": "time (UTC)"}
    )

    ds["sample"].attrs.update({"units": "1", "long_name": "Sample in burst"})

    ds["Tx_1211"].attrs.update(
        {
            "units": "degree_C",
            "long_name": "Instrument Internal Temperature",
            # "standard_name": "sea_water_temperature",
            "epic_code": "1211",
        }
    )

    if "ea" in ds.attrs["instrument_type"]:
        ds.attrs["instrument_type"] = "EofE ECHOLOGGER EA400 profiling altimeter"

        ds["AGC_1202"].attrs.update(
            {
                "units": "counts",
                "long_name": "Average Echo Intensity",
                # "generic_name": "AGC",
                "epic_code": "1202",
            }
        )

        ds["Ptch_1216"].attrs.update(
            {
                "units": "degrees",
                "long_name": "Instrument Pitch",
                "standard_name": "platform_pitch",
                "epic_code": "1216",
            }
        )

        ds["Roll_1217"].attrs.update(
            {
                "units": "degrees",
                "long_name": "Instrument Roll",
                "standard_name": "platform_roll",
                "epic_code": "1217",
            }
        )

    if "aa" in ds.attrs["instrument_type"]:
        ds.attrs["instrument_type"] = "EofE ECHOLOGGER AA400 altimeter"

        ds["AMP_723"].attrs.update(
            {
                "units": "percent",
                "long_name": "Acoustic Signal Amplitude Strength",
                "epic_code": "723",
            }
        )

    return ds


def calc_bin_height(ds):
    # modified from qaqc.check_orientation

    print("Calculating center of bin distance from transducer")

    # np.linspace(start,stop,num)
    # start: 0, because first bin is at transducer (confirmed with EofE), add (ds.attrs["Bin_size_m"] / 2) for center of bin as point of reference
    # stop: number of bins - 1 * bin size, add (ds.attrs["Bin_size_m"] / 2) for center of bin as point of reference
    # num: number of bins (Bin_count)

    prefix = ds.attrs["instrument_type"].upper()

    ds["bindist"] = xr.DataArray(
        np.linspace(
            0 + (ds.attrs[f"{prefix}Bin_size_m"] / 2),
            (
                ((ds.attrs[f"{prefix}Bin_count"] - 1) * ds.attrs[f"{prefix}Bin_size_m"])
                + (ds.attrs[f"{prefix}Bin_size_m"] / 2)
            ),
            num=ds.attrs[f"{prefix}Bin_count"],
        ),
        dims="bins",
    )

    print(
        "Calculating center of bin height from seafloor as: initial instrument height +/- bin(center) distance from transducer"
    )

    if ds.attrs["orientation"].upper() == "DOWN":
        ds["bin_height"] = (
            ds.attrs["initial_instrument_height"] - ds["bindist"]
        )  # get bin distance referenced from sea floor

        math_sign = "-"

    elif ds.attrs["orientation"].upper() == "UP":
        ds["bin_height"] = (
            ds.attrs["initial_instrument_height"] + ds["bindist"]
        )  # get bin distance referenced from sea floor

        math_sign = "+"

    # add attributes
    ds["bindist"].attrs.update(
        {
            "units": "m",
            "long_name": "bin(center) distance from transducer",
            "positive": "%s" % ds.attrs["orientation"],
            "note": "Distance is along profile from instrument head to center of bin",
        }
    )

    ds["bin_height"].attrs.update(
        {
            "units": "m",
            "long_name": "bin(center) distance from seafloor",
            "positive": "up",
            "note": f"Distance is along profile from seafloor to center of bin. Calculated as initial instrument height {math_sign} bin(center) distance from transducer based on {ds.attrs[f'{prefix}SoundSpeed_mps']} m/s sound vel.",
            # % math_sign,
        }
    )

    # round to mm
    ds["bindist"] = ds["bindist"].round(3)
    ds["bin_height"] = ds["bin_height"].round(3)

    return ds


def calc_cor_brange(ds):
    print("Correcting distance to bed (brange) using adjusted sound speed")
    # here the brange is still called Altitude_m but variable name will change to brange later in code

    prefix = ds.attrs["instrument_type"].upper()

    # Correct using adjusted sound speed
    # distance (m) = time (sec) x sound speed (m/sec)
    time_sec = xr.DataArray((ds.Altitude_m) / ds.attrs[f"{prefix}SoundSpeed_mps"])
    #
    # seawater.svel(s,t,p); s = average salinity (psu) from exo, t = temp (c) from altimeter, p = approximate pressure (db) calculated from instrument depth +/- median(Altitude)/2
    if ds.attrs["orientation"].upper() == "DOWN":
        p = (
            ds.attrs["WATER_DEPTH"]
            - ds.attrs["initial_instrument_height"]
            + ds["Altitude_m"].median() / 2
        )
        math_sign = "-"
    elif ds.attrs["orientation"].upper() == "UP":
        p = (
            ds.attrs["WATER_DEPTH"]
            - ds.attrs["initial_instrument_height"]
            - ds["Altitude_m"].median() / 2
        )
        math_sign = "+"

    soundspd = gsw.sound_speed(ds.attrs["average_salinity"], ds.Temperature_C, p)
    ds["brange"] = xr.DataArray(time_sec * soundspd).round(3)  # round brange to mm

    histtext = (
        f"Adjusted sound velocity calculated using sound_speed(s,t,p) from gsw toolbox (https://teos-10.github.io/GSW-Python/). Inputs: Salinity (s) from average salinity of "
        f"{ds.attrs['average_salinity']} PSU, temperature (t) from ea400 internal temperature measurements, pressure (p) from instrument depth {math_sign} median(altitude)/2. "
    )

    ds = utils.insert_history(ds, histtext)

    ds["brange"].attrs.update(
        {
            "units": "m",
            "long_name": "Altimeter range to boundary",
            "standard_name": "altimeter_range",
            "note": "Calculated using adjusted speed of sound",
        }
    )

    return ds


def calc_boundary_elev(ds):
    # find seabed elevation referenced to datum

    if "NAVD88_ref" in ds.attrs:
        ds.attrs["geopotential_datum_name"] = "NAVD88"

        print(
            "Calculating seabed elevation on %s datum"
            % ds.attrs["geopotential_datum_name"]
        )

        if ds.attrs["orientation"].upper() == "DOWN":
            ds["boundary_elevation"] = xr.DataArray(
                ds.attrs["NAVD88_ref"]
                + (ds.brange * -1)
                + ds.attrs["initial_instrument_height"]
            )

        elif ds.attrs["orientation"].upper() == "UP":
            ds["boundary_elevation"] = xr.DataArray(
                ds.attrs["NAVD88_ref"]
                + ds.brange
                + ds.attrs["initial_instrument_height"]
            )

    elif "height_above_geopotential_datum" in ds.attrs:
        print(
            "Calculating seabed elevation on %s datum"
            % ds.attrs["geopotential_datum_name"]
        )

        if ds.attrs["orientation"].upper() == "DOWN":
            ds["boundary_elevation"] = xr.DataArray(
                ds.attrs["height_above_geopotential_datum"]
                + (ds.brange * -1)
                + ds.attrs["initial_instrument_height"]
            )

        elif ds.attrs["orientation"].upper() == "UP":
            ds["boundary_elevation"] = xr.DataArray(
                ds.attrs["height_above_geopotential_datum"]
                + ds.brange
                + ds.attrs["initial_instrument_height"]
            )

    else:
        ds.attrs["geopotential_datum_name"] = "LMSL"
        print(
            "Calculating seabed elevation on %s datum"
            % ds.attrs["geopotential_datum_name"]
        )

        if ds.attrs["orientation"].upper() == "DOWN":
            ds["boundary_elevation"] = xr.DataArray(
                ds.attrs["WATER_DEPTH"]
                + ds.brange
                - ds.attrs["initial_instrument_height"]
            )

        if ds.attrs["orientation"].upper() == "UP":
            ds["boundary_elevation"] = xr.DataArray(
                ds.attrs["WATER_DEPTH"]
                + (ds.brange * -1)
                - ds.attrs["initial_instrument_height"]
            )

    if "height_above_geopotential_datum" in ds.attrs or "NAVD88_ref" in ds.attrs:
        ds["boundary_elevation"].attrs.update(
            {
                "units": "m",
                "long_name": f"boundary height referenced to {ds.attrs['geopotential_datum_name']} datum",
                # % ds.attrs["geopotential_datum_name"],
                "standard_name": "height_above_geopotential_datum",
                "positive": "up",
                "note": f"Corrected brange from adjusted sound speed and referenced to {ds.attrs['geopotential_datum_name']} datum",
                # % ds.attrs["geopotential_datum_name"],
            }
        )

    elif ds.attrs["geopotential_datum_name"] == "LMSL":
        ds["boundary_elevation"].attrs.update(
            {
                "units": "m",
                "long_name": f"boundary depth referenced to {ds.attrs['geopotential_datum_name']} datum",
                # % ds.attrs["geopotential_datum_name"],
                "standard_name": "depth_below_geoid",
                "positive": "down",
                "note": f"Corrected brange from adjusted sound speed and referenced to {ds.attrs['geopotential_datum_name']} datum",
                # % ds.attrs["geopotential_datum_name"],
            }
        )

    ds["boundary_elevation"] = ds["boundary_elevation"].round(
        3
    )  # round seabed_elevation to mm

    histtext = (
        "add boundary_elevation using speed of sound corrected brange and ref datum"
    )
    ds = utils.insert_history(ds, histtext)

    return ds


def calc_cor_bin_height(ds):
    print("Calculating corrected bin height from adjusted speed of sound")

    prefix = ds.attrs["instrument_type"].upper()

    binheight = (
        xr.DataArray(ds.bin_height)
        .expand_dims({"time": ds.time, "sample": ds.sample})
        .transpose("bins", "time", "sample")
    )
    # calculate travel time per bin
    time_sec = xr.DataArray((binheight) / ds.attrs[f"{prefix}SoundSpeed_mps"])
    # calculate sound speed

    p = (
        ds.attrs["WATER_DEPTH"] - binheight
    )  # use water depth and bin_height to estimate p

    # spd = sw.svel(ds.attrs["average_salinity"], ds.Temperature_C, p)
    # speed up finding sound speed by taking mean across sample dim then expand dims after
    spd2 = gsw.sound_speed(
        ds.attrs["average_salinity"],
        ds["Temperature_C"].mean(dim="sample").values,
        p.mean(dim="sample").values,
    )
    soundspd = (
        xr.DataArray(spd2, dims=["bins", "time"])
        .expand_dims({"sample": ds.sample})
        .transpose("bins", "time", "sample")
    )

    ds["cor_bin_height"] = xr.DataArray(time_sec * soundspd).transpose(
        "time", "sample", "bins"
    )

    ds["cor_bin_height"].attrs.update(
        {
            "units": "m",
            "long_name": "corrected bin(center) distance from seafloor",
            "positive": "up",
            "note": "Distance is along profile from seafloor to center of bin based on adjusted sound vel.",
        }
    )

    # round to mm
    ds["cor_bin_height"] = ds["cor_bin_height"].round(3)

    return ds


def average_burst(ds):
    ds = ds.mean(
        "sample", skipna=True, keep_attrs=True
    )  # take mean across 'sample' dim

    # round brange and seabed_elevation to 3 decimal places (mm)
    if "brange" in ds:
        ds["brange"] = ds["brange"].round(3)

    if "seabed_elevation" in ds:
        ds["seabed_elevation"] = ds["seabed_elevation"].round(3)

    return ds


def read_aa_instmet(basefile):
    with open(basefile) as f:
        instmeta = {}
        row = ""
        while "   Date       Time" not in row:
            row = f.readline().rstrip()
            dat = row.split()
            if " Device" in row:
                instmeta["serial_number"] = dat[2]
            elif "- range" in row:
                instmeta["AARange_m"] = float(dat[3]) / 1000
            elif "- interval" in row:
                instmeta["AAInterval_between_series_sec"] = float(dat[3])
            elif "- series" in row:
                instmeta["AAPulses_in_series_num"] = int(dat[3])
            elif "- threshold" in row:
                instmeta["AAThreshold_percent"] = int(dat[3])
            elif "- offset" in row:
                instmeta["AAOffset_m"] = float(dat[3]) / 1000
            elif "- deadzone" in row:
                instmeta["AADeadzone_m"] = float(dat[3]) / 1000
            elif "- txlength" in row:
                instmeta["AAPulseLength_microsec"] = int(dat[3])
            elif "- amplitude" in row:
                instmeta["AAAmplitude_Threshold_percent"] = int(dat[3])
            elif "-" and "sampling" in row:
                instmeta["AASampling_rate_Hz"] = int(dat[3])
            elif "Total Number of Series" in row:
                instmeta["AANSeries"] = int(dat[5])
            elif "Total Number of Records" in row:
                instmeta["AANRecords"] = int(dat[5])

    instmeta["AASoundSpeed_mps"] = 1500  # aa400 fixed speed of sound value
    return instmeta


def load_aa_point(filnam, instmeta, skiprows=None, skipfooter=None, encoding="utf-8"):
    """Read data from an EofE AA400 .log file into an xarray
    Dataset.

    Parameters
    ----------
    filnam : string
        The filename
    skiprows : int, optional
        How many header rows to skip. Default None
    skipfooter : int, optional
        How many footer rows to skip. Default None
    encoding : string, optional
        File encoding. Default 'utf-8'
    Returns
    -------
    xarray.Dataset
        An xarray Dataset of the altimeter data
    """

    df = pd.read_csv(
        filnam,
        skiprows=skiprows,
        skipfooter=skipfooter,
        header=None,
        usecols=np.arange(0, 7),  # Use first 7 columns of data (altitude in mm)
        sep=r"[, ]+",  # Delimiter is comma or whitespace
        names=[
            "Date",
            "Time",
            "Ping",
            "Altitude_m",
            "Temperature_C",
            "Battery_mV",
            "AmplitudeFS",
        ],
        encoding=encoding,
        engine="python",
        index_col=False,
    )

    # Need to use because parse_dates in pd.read_csv is deprecated
    df["time"] = pd.to_datetime(
        df["Date"].astype(str) + (df["Time"]), format="%Y%m%d%H:%M:%S.%f"
    )
    df.pop("Date")
    df.pop("Time")

    # Convert altitude from mm to m
    df["Altitude_m"] = df["Altitude_m"] / 1000

    # Convert to dictionary and reshape
    point = df.to_dict(orient="list")
    samples = instmeta["AAPulses_in_series_num"]

    for k in point:
        point[k] = np.reshape(point[k], (-1, samples))

    time = point["time"][:, 0]

    ds = xr.Dataset()
    ds["time"] = xr.DataArray(pd.to_datetime(time), dims="time")
    ds["sample"] = xr.DataArray(np.arange(0, samples), dims="sample")

    # add variables to xarray
    for k in point:
        if "time" not in k:
            ds[k] = xr.DataArray(point[k], dims=("time", "sample"))

    # need to retype to float
    ds["AmplitudeFS"] = ds["AmplitudeFS"].astype(dtype="float")

    return ds


def trim_alt(ds, data_vars=["Altitude_m", "Counts", "Temperature_C"]):
    """Trim altimeter data when out of water during a deployment and by bin range if specified

    Parameters
    ----------
    ds : xarray.Dataset
        The xarray Dataset
    data_vars : array_like
        List of variables to trim. Default ['Altitude_m', 'Counts', 'Temperature_C'].

    Returns
    -------
    xarray.Dataset
        Dataset with trimmed data
    """

    prefix = ds.attrs["instrument_type"].upper()

    if "trim_method" in ds.attrs:
        trm_list = ds.attrs["trim_method"]

        if not isinstance(trm_list, list):  # make sure it is a list before looping
            trm_list = [trm_list]

        for trm_meth in trm_list:
            if trm_meth.lower() == "altitude":
                print("Trimming using altitude data")
                # need to use altitude values before starting trimming
                for var in data_vars:
                    ds[var] = ds[var].where(
                        ds["Altitude_m"] >= ds.attrs[f"{prefix}Deadzone_m"]
                    )
                    ds[var] = ds[var].where(
                        ds["Altitude_m"] <= ds.attrs[f"{prefix}Range_m"]
                    )
                    print(f"Trimming {var}")

                histtext = "Trimmed altimeter data using Altimeter_m = 0."

                ds = utils.insert_history(ds, histtext)

            elif trm_meth.lower() == "bin range":
                print("Trimming using good_bins of %s" % str(ds.attrs["good_bins"]))
                if "bins" in ds.coords:
                    # trim coordinate bins
                    ds = ds.isel(
                        bins=slice(ds.attrs["good_bins"][0], ds.attrs["good_bins"][1])
                    )
                    # reset Bin_count attribute
                    ds.attrs["Bin_count"] = (
                        ds.attrs["good_bins"][1] - ds.attrs["good_bins"][0]
                    )

                histtext = (
                    "Removed extra bins from altimeter data using good_bins attribute."
                )

                ds = utils.insert_history(ds, histtext)

            else:
                print("Did not trim altimeter data")

    return ds
