import numpy as np
import pandas as pd
import xarray as xr
import seawater as sw

from .core import utils


def log_to_cdf(metadata):

    basefile = metadata["basefile"]

    if "prefix" in metadata:
        basefile = metadata["prefix"] + basefile

    # utils.check_valid_metadata(metadata)

    # get instrument metadata from the LOG file
    # Assume ea400 in burst mode if no attrs
    if "instrument_type" not in metadata or metadata["instrument_type"] == "ea":
        if "instrument_type" not in metadata:
            print("instrument_type not specified assumed to be type == ea")
            metadata["instrument_type"] = "ea"

        instmeta = read_ea_instmet(basefile + ".log")
        metadata["instmeta"] = instmeta
        print("Loading LOG file")
        # load point sensor data (temp, altitude, pitch, roll, ping #, sample #)
        ds = load_ea_point(basefile + ".log", metadata)

        ds = utils.write_metadata(ds, metadata)

        # load profile sensor data (counts)
        ds = load_ea_profile(ds, basefile + ".log")

    # Else, check for aa400
    elif metadata["instrument_type"] == "aa":
        instmeta = read_aa_instmet(basefile + ".log")
        metadata["instmeta"] = instmeta
        print("Loading LOG file")
        # load point sensor data (temp, altitude, ping #)
        ds = load_aa_point(basefile + ".log", metadata)
        ds = utils.write_metadata(ds, metadata)
    else:
        raise ValueError(
            "instrument_type in config file, {:s}, is invalid".format(
                metadata["instrument_type"]
            )
        )

    ds = utils.shift_time(ds, 0)

    # configure file
    print("Configuring .cdf file")

    if "prefix" in ds.attrs:
        cdf_filename = ds.attrs["prefix"] + ds.attrs["filename"] + "-raw.cdf"
    else:
        cdf_filename = ds.attrs["filename"] + "-raw.cdf"

    ds["time"] = ds["time"].astype("datetime64[s]")
    ds.to_netcdf(
        cdf_filename, unlimited_dims=["time"], encoding={"time": {"dtype": "i4"}}
    )

    print("Finished writing data to %s" % cdf_filename)

    return ds


def cdf_to_nc(cdf_filename):
    """
    Load a "raw" .cdf file and generate a processed .nc file
    """

    # Load raw .cdf data
    ds = xr.load_dataset(cdf_filename)

    # create burst num variable
    ds = burst_num(ds)

    # Clip data to in/out water times or via good_ens
    ds = utils.clip_ds(ds)

    # calculate bin height for profiling echologger (i.e. type == ea)
    if "bins" in ds:
        calc_bin_height(ds)

    # calculate corrected altitude (distance to bed/b_range) with adjusted sound speed
    ds = calc_cor_brange(ds)  # 12/23/21

    # calculate corrected bin height (on NAVD88 datum) with adjusted sound speed
    if "bins" in ds:
        ds = calc_cor_bin_height(ds)

    ds = calc_seabed_elev(ds)

    # if "bins" in ds:
    #   ds = utils.create_z_bindist(ds)  #create z for profile
    # else:
    ds = utils.create_z(ds)

    # swap bin dim with bin_height
    ds = ds_swap_dims(ds)  # swap vert dim to z

    # rename variables
    ds = ds_rename_vars(ds)

    # drop some vars after renaming
    for k in [
        "bins",
        "ping",
        "ping_num_in_series",
        "Altitude_m",
        "Battery_mV",
    ]:
        if k in ds:
            ds = ds.drop_vars(k)

    # add lat/lons as coordinates
    ds = utils.ds_add_lat_lon(ds)

    # add attributes to each variable
    ds = ds_add_attrs(ds)

    # add metadata to global atts
    ds = utils.add_start_stop_time(ds)
    ds = utils.add_delta_t(ds)

    # Write to .nc file
    print("Writing cleaned/trimmed burst data and averaged burst data to .nc file")
    nc_filename = ds.attrs["filename"] + "b-cal.nc"

    ds.to_netcdf(
        nc_filename, unlimited_dims=["time"], encoding={"time": {"dtype": "i4"}}
    )
    print("Done writing netCDF file", nc_filename)

    # Average busrt and write to -a.nc file 12/23/21
    ds = average_burst(ds)

    # assign min/max
    ds = utils.add_min_max(ds)

    nc_filename = ds.attrs["filename"] + "-a.nc"

    # ds['time']=ds['time'].astype('datetime64[s]')
    ds.to_netcdf(
        nc_filename, unlimited_dims=["time"], encoding={"time": {"dtype": "i4"}}
    )  # 5/20/22- save time as int32

    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])

    print("Done writing burst averaged netCDF file", nc_filename)

    return ds


def read_ea_instmet(basefile):
    with open(basefile, "r") as f:

        instmeta = {}
        row = ""
        while "##DataStart" not in row:
            row = f.readline().rstrip()
            dat = row.split()
            if "#DeviceID" in row:
                instmeta["DeviceID"] = row[10:]
                instmeta["serial_number"] = dat[1]
            elif "#NSamples" in row:
                instmeta["Bin_count"] = int(dat[1])
            elif "#Resolution,m" in row:
                instmeta["Bin_size_m"] = float(dat[1])
            elif "#SoundSpeed,mps" in row:
                instmeta["SoundSpeed_mps"] = float(dat[1])
            elif "#Tx_Frequency,Hz" in row:
                instmeta["Tx_Frequency_Hz"] = float(dat[1])
            elif "#Range,m" in row:
                instmeta["Range_m"] = float(dat[1])
            elif "#Pulse period,sec" in row:
                instmeta["Pulse_period_sec"] = float(dat[2])
            elif "#Pulses in series,num" in row:
                instmeta["Pulses_in_series_num"] = int(dat[3])
            elif "#Interval between series,sec" in row:
                instmeta["Interval_between_series_sec"] = float(dat[3])
            elif "#Threshold,%" in row:
                instmeta["Threshold_percent"] = int(dat[1])
            elif "#Offset,m" in row:
                instmeta["Offset_m"] = float(dat[1])
            elif "#Deadzone,m" in row:
                instmeta["Deadzone_m"] = float(dat[1])
            elif "#PulseLength,uks" in row:
                instmeta["PulseLength_microsec"] = float(dat[1])
            elif "#TVG_Gain,dB" in row:
                instmeta["TVG_Gain_dB"] = float(dat[1])
            elif "#TVG_Slope,dB/km" in row:
                instmeta["TVG_Slope_dBkm"] = float(dat[1])
            elif "#TVG_Mode" in row:
                instmeta["TVG_Mode"] = int(dat[1])
            elif "#OutputMode" in row:
                instmeta["OutputMode"] = int(dat[1])

    return instmeta


def load_ea_point(basefile, metadata):
    with open(basefile, "r") as f:

        data = f.read().splitlines()

        point = {}

        TimeLocal = []
        TimeUTC = []
        Ping = []
        Ping_num_in_series = []
        Altitude_m = []
        Temperature_C = []
        Roll_deg = []
        Pitch_deg = []

        for row in data:
            dat = row.split()
            if "#TimeLocal" in row:
                TimeLocal.append(dat[1] + " " + dat[2])
            elif "#TimeUTC" in row:
                TimeUTC.append(dat[1] + " " + dat[2])
            elif "#Ping  " in row:
                Ping.append(float(dat[1]))
            elif "#Ping num in series" in row:
                Ping_num_in_series.append(float(dat[4]))
            elif "#Altitude,m" in row:
                Altitude_m.append(float(dat[1]))
            elif "#Temperature" in row:
                Temperature_C.append(float(dat[1]))
            elif "#Pitch,deg" in row:
                Pitch_deg.append(float(dat[1]))
            elif "#Roll,deg" in row:
                Roll_deg.append(float(dat[1]))

    # Add lists to point dictionary
    point["TimeLocal"] = TimeLocal
    point["TimeUTC"] = TimeUTC
    point["Ping"] = Ping
    point["Ping_num_in_series"] = Ping_num_in_series
    point["Altitude_m"] = Altitude_m
    point["Temperature_C"] = Temperature_C
    point["Pitch_deg"] = Pitch_deg
    point["Roll_deg"] = Roll_deg

    for k in point:
        point[k] = np.array(point[k])

    point["TimeUTC"] = pd.to_datetime(point["TimeUTC"])
    point["TimeLocal"] = pd.to_datetime(point["TimeLocal"])

    # reshape point data
    samples = metadata["instmeta"]["Pulses_in_series_num"]
    n = metadata["instmeta"]["Bin_count"]
    for k in point:
        if "Time" not in k:
            point[k] = point[k].reshape((-1, samples))

    time = point["TimeUTC"][::samples]
    ds = xr.Dataset()
    ds["time"] = xr.DataArray(time, dims="time")
    ds["sample"] = xr.DataArray(np.arange(0, samples), dims="sample")
    ds["bins"] = xr.DataArray(np.arange(n), dims="bins")

    # add dimensions to variables
    for k in point:
        if "Time" not in k:
            ds[k] = xr.DataArray(point[k], dims=("time", "sample"))

    return ds


def load_ea_profile(ds, basefile):

    profile = []

    with open(basefile) as f:  # read in profile echo data
        for row in f:
            if row.rstrip() == "##DataStart":
                for row in f:
                    if row.rstrip() == "##DataEnd":
                        break
                    profile.append(float(row))

    # reshape profile data
    samples = ds.Pulses_in_series_num
    n = ds.Bin_count
    profile = np.array(profile)
    profile = profile.reshape((-1, samples, n))
    ds["Counts"] = xr.DataArray(
        profile,
        dims=("time", "sample", "bins"),
        coords=[ds["time"], ds["sample"], ds["bins"]],
    )

    return ds


def burst_num(ds):

    ds["burst"] = xr.DataArray(
        np.arange(1, len(ds["time"]) + 1, 1, dtype="int32"), dims="time"
    )

    return ds


def ds_rename_vars(ds):
    # modified from exo.ds_rename_vars
    varnames = {
        "Ping": "ping",
        "Ping_num_in_series": "ping_num_in_series",
        "Temperature_C": "Tx_1211",  # 12/23/21
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


def ds_add_attrs(ds):  # 12/23/21
    # modified from exo.ds_add_attrs
    ds = utils.ds_coord_no_fillvalue(ds)

    ds["time"].attrs.update({"standard_name": "time", "axis": "T"})

    ds["sample"].attrs.update({"units": "1", "long_name": "Sample in burst"})

    ds["burst"].attrs.update(
        {
            "units": "1",
            "long_name": "Burst number",
            "generic_name": "record",
            "epic_code": "1207",
            "coverage_content_type": "physicalMeasurement",
        }
    )

    ds["Tx_1211"].attrs.update(
        {
            "units": "degree_C",
            "long_name": "Instrument Internal Temperature",
            "standard_name": "sea_water_temperature",
            "epic_code": "1211",
        }
    )

    if "ea" in ds.attrs["instrument_type"]:
        ds["AGC_1202"].attrs.update(
            {
                "units": "counts",
                "long_name": "Average Echo Intensity",
                "generic_name": "AGC",
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
        ds["AMP_723"].attrs.update(
            {
                "units": "percent",
                "long_name": "Acoustic Signal Amplitude Strength",
                "epic_code": "723",
            }
        )

    # add initial height information and fill values to variabels
    def add_attributes(var, dsattrs):
        if "ea" in dsattrs["instrument_type"]:
            var.attrs.update(
                {
                    "initial_instrument_height": dsattrs["initial_instrument_height"],
                    "height_depth_units": "m",
                    "sensor_type": "ECHOLOGGER EA400",
                }
            )
        elif "aa" in dsattrs["instrument_type"]:
            var.attrs.update(
                {
                    "initial_instrument_height": dsattrs["initial_instrument_height"],
                    "height_depth_units": "m",
                    "sensor_type": "ECHOLOGGER AA400",
                }
            )

    # don't include all attributes for coordinates that are also variables
    for var in ds.variables:
        if (var not in ds.coords) and ("time" not in var):
            add_attributes(ds[var], ds.attrs)

    return ds


def calc_bin_height(ds):
    # modified from qaqc.check_orientation

    print("Calculating center of bin distance from transducer")

    # np.linspace(start,stop,num)
    ##start: 0, because first bin is at transducer (confirmed with EofE), add (ds.attrs["Bin_size_m"] / 2) for center of bin as point of reference
    ##stop: number of bins - 1 * bin size, add (ds.attrs["Bin_size_m"] / 2) for center of bin as point of reference
    ##num: number of bins (Bin_count)

    ds["bindist"] = xr.DataArray(
        np.linspace(
            0 + (ds.attrs["Bin_size_m"] / 2),  # changed ds.attrs["Bin_size_m"] to 0
            (
                (
                    (ds.attrs["Bin_count"] - 1) * ds.attrs["Bin_size_m"]
                )  # deleted ds.attrs["Bin_size_m"] +
                + (ds.attrs["Bin_size_m"] / 2)
            ),
            num=ds.attrs["Bin_count"],
        ),
        dims="bins",
    )

    print(
        "Calculating center of bin height from seafloor as: initial intrument height - bin(center) distance from transducer"
    )

    if ds.attrs["orientation"] == "down" or ds.attrs["orientation"] == "DOWN":

        ds["bin_height"] = (
            ds.attrs["initial_instrument_height"] - ds["bindist"]
        )  # get bin distance referenced from sea floor

        math_sign = "-"

    elif ds.attrs["orientation"] == "up" or ds.attrs["orientation"] == "UP":

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
            "note": "Distance is along profile from seafloor to center of bin. Calculated as initial instrument height %s bin(center) distance from transducer based on 1500 m/s sound vel."
            % math_sign,
        }
    )

    # round to mm
    ds["bindist"] = ds["bindist"].round(3)
    ds["bin_height"] = ds["bin_height"].round(3)

    return ds


def calc_cor_brange(ds):
    print("Correcting distance to bed (brange) using adjusted sound speed")
    # here the brange is still called Altitude_m but variable name will change to brange later in code

    # Correct using adjusted sound speed
    # distance (m) = time (sec) x sound speed (m/sec)
    time_sec = xr.DataArray((ds.Altitude_m) / ds.attrs["SoundSpeed_mps"])

    # seawater.svel(s,t,p); s = average salinity (psu) from exo, t = temp (c) from altimeter, p = approximate pressure (db) calculated from Altitude_m
    soundspd = sw.svel(ds.attrs["average_salinity"], ds.Temperature_C, ds.Altitude_m)
    ds["brange"] = xr.DataArray(time_sec * soundspd).round(3)  # round brange to mm

    histtext = "Adjusted sound velocity calculated using svel(s,t,p) from seawater toolbox (https://pythonhosted.org/seawater/eos80.html#seawater.eos80.svel). Svel inputs: Salinity (s) from average salinity of %s PSU, temperature (t) from ea400 internal temperature measurements, pressure (p) from raw ea400 altitude measurements. "

    ds = utils.insert_history(ds, histtext)

    ds["brange"].attrs.update(
        {
            "units": "m",
            "long_name": "sensor range to boundary",
            "standard_name": "height_above_sea_floor",
            "note": "Calculated using adjusted speed of sound",
        }
    )

    return ds


def calc_seabed_elev(ds):
    # find seabed elevation referenced to datum
    histtext = (
        "add seabed_elevation using speed of sound corrected brange and ref datum"
    )
    ds = utils.insert_history(ds, histtext)

    if "NAVD88_ref" in ds.attrs:
        ds.attrs["geopotential_datum_name"] = "NAVD88"

        print(
            "Calculating seabed elevation on %s datum"
            % ds.attrs["geopotential_datum_name"]
        )

        if ds.attrs["orientation"] == "DOWN" or ds.attrs["orientation"] == "down":

            ds["seabed_elevation"] = xr.DataArray(
                ds.attrs["NAVD88_ref"]
                + (ds.brange * -1)
                + ds.attrs["initial_instrument_height"]
            )

        elif ds.attrs["orientation"] == "UP" or ds.attrs["orientation"] == "up":

            ds["seabed_elevation"] = xr.DataArray(
                ds.attrs["NAVD88_ref"]
                + ds.brange
                + ds.attrs["initial_instrument_height"]
            )

    elif "height_above_geopotential_datum" in ds.attrs:
        print(
            "Calculating seabed elevation on %s datum"
            % ds.attrs["geopotential_datum_name"]
        )

        if ds.attrs["orientation"] == "DOWN" or ds.attrs["orientation"] == "down":

            ds["seabed_elevation"] = xr.DataArray(
                ds.attrs["height_above_geopotential_datum"]
                + (ds.brange * -1)
                + ds.attrs["initial_instrument_height"]
            )

        elif ds.attrs["orientation"] == "UP" or ds.attrs["orientation"] == "up":

            ds["seabed_elevation"] = xr.DataArray(
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

        if ds.attrs["orientation"] == "DOWN" or ds.attrs["orientation"] == "down":

            ds["seabed_elevation"] = xr.DataArray(
                ds.attrs["WATER_DEPTH"]
                + ds.brange
                - ds.attrs["initial_instrument_height"]
            )

        if ds.attrs["orientation"] == "UP" or ds.attrs["orientation"] == "up":

            ds["seabed_elevation"] = xr.DataArray(
                ds.attrs["WATER_DEPTH"]
                + (ds.brange * -1)
                + ds.attrs["initial_instrument_height"]
            )

    if "height_above_geopotential_datum" in ds.attrs or "NAVD88_ref" in ds.attrs:
        ds["seabed_elevation"].attrs.update(
            {
                "units": "m",
                "long_name": "seafloor height referenced to %s datum"
                % ds.attrs["geopotential_datum_name"],
                "standard_name": "height_above_geopotential_datum",
                "positive": "up",
                "note": "Corrected brange from adjusted sound speed and referenced to %s datum"
                % ds.attrs["geopotential_datum_name"],
            }
        )

    elif ds.attrs["geopotential_datum_name"] == "LMSL":
        ds["seabed_elevation"].attrs.update(
            {
                "units": "m",
                "long_name": "sea floor depth referenced to %s datum"
                % ds.attrs["geopotential_datum_name"],
                "standard_name": "sea_floor_depth_below_mean_sea_level",
                "positive": "down",
                "note": "Corrected brange from adjusted sound speed and referenced to %s datum"
                % ds.attrs["geopotential_datum_name"],
            }
        )

    ds["seabed_elevation"] = ds["seabed_elevation"].round(
        3
    )  # round seabed_elevation to mm

    return ds


def calc_cor_bin_height(ds):
    print("Calculating corrected bin height from adjusted speed of sound")
    binheight = (
        xr.DataArray(ds.bin_height)
        .expand_dims({"time": ds.time, "sample": ds.sample})
        .transpose("bins", "time", "sample")
    )
    # calculate travel time per bin
    time_sec = xr.DataArray((binheight) / ds.attrs["SoundSpeed_mps"])
    # calculate sound speed
    spd = sw.svel(ds.attrs["average_salinity"], ds.Temperature_C, binheight)
    ds["cor_bin_height"] = xr.DataArray(time_sec * spd).transpose(
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
    ds["burst"] = ds.burst.astype(
        dtype="int32"
    )  # need to retype to int32 bc np/xarray changes int to float when averaging

    # round brange and seabed_elevation to 3 decimal places (mm)
    if "brange" in ds:
        ds["brange"] = ds["brange"].round(3)

    if "seabed_elevation" in ds:
        ds["seabed_elevation"] = ds["seabed_elevation"].round(3)

    return ds


def ds_swap_dims(ds):  # swap vert dim to z
    # need to preserve z attrs because swap_dims will remove them
    attrsbak = ds["z"].attrs
    for v in ds.data_vars:
        if "bins" in ds[v].coords:
            ds[v] = ds[v].swap_dims({"bins": "z"})

    ds["z"].attrs = attrsbak

    return ds


def read_aa_instmet(basefile):
    with open(basefile, "r") as f:

        instmeta = {}
        row = ""
        while "   Date       Time" not in row:
            row = f.readline().rstrip()
            dat = row.split()
            if " Device" in row:
                instmeta["DeviceID"] = dat[2]
            elif "- range" in row:
                instmeta["Range_m"] = float(dat[3]) / 1000
            elif "- interval" in row:
                instmeta["Interval_between_series_sec"] = float(dat[3])
            elif "- series" in row:
                instmeta["Pulses_in_series_num"] = int(dat[3])
            elif "- threshold" in row:
                instmeta["Threshold_percent"] = int(dat[3])
            elif "- offset" in row:
                instmeta["Offset_m"] = float(dat[3]) / 1000
            elif "- deadzone" in row:
                instmeta["Deadzone_m"] = float(dat[3]) / 1000
            elif "- txlength" in row:
                instmeta["PulseLength_microsec"] = int(dat[3])
            elif "- amplitude" in row:
                instmeta["Amplitude_Threshold_percent"] = int(dat[3])
            elif "-" and "sampling" in row:
                instmeta["Sampling_rate_Hz"] = int(dat[3])
            elif "Total Number of Series" in row:
                instmeta["NSeries"] = int(dat[5])
            elif "Total Number of Records" in row:
                instmeta["NRecords"] = int(dat[5])

    instmeta["SoundSpeed_mps"] = 1500  # aa400 fixed speed of sound value
    return instmeta


def load_aa_point(basefile, metadata):
    with open(basefile, "r") as f:

        if "skiprows" in metadata:
            for k in np.arange(0, metadata["skiprows"]):
                line = f.readline()

        else:
            line = ""
            while "   Date       Time" not in line:
                line = f.readline()
            line = f.readline()

        data = f.read().splitlines()

        point = {}

        TimeUTC = []
        Ping = []
        Altitude_m = []
        Temperature_C = []
        Battery_mV = []
        AmplitudeFS = []

        for row in data:
            dat = np.array(row.split())
            TimeUTC.append(dat[0] + " " + dat[1])
            Ping.append(float(dat[2]))
            Altitude_m.append(float(dat[3]) / 1000)
            Temperature_C.append(float(dat[4]))
            Battery_mV.append(float(dat[5]))
            AmplitudeFS.append(float(dat[6]))

    # Add lists to point dictionary
    point["TimeUTC"] = TimeUTC
    point["Ping"] = Ping
    point["Altitude_m"] = Altitude_m
    point["Temperature_C"] = Temperature_C
    point["Battery_mV"] = Battery_mV
    point["AmplitudeFS"] = AmplitudeFS

    for k in point:
        point[k] = np.array(point[k])

    point["TimeUTC"] = pd.to_datetime(point["TimeUTC"])

    # reshape point data
    samples = metadata["instmeta"]["Pulses_in_series_num"]

    for k in point:
        if "Time" not in k:
            point[k] = point[k].reshape((-1, samples))

    time = point["TimeUTC"][::samples]
    ds = xr.Dataset()
    ds["time"] = xr.DataArray(time, dims="time")
    ds["sample"] = xr.DataArray(np.arange(0, samples), dims="sample")
    # add dimensions to variables
    for k in point:
        if "Time" not in k:
            ds[k] = xr.DataArray(point[k], dims=("time", "sample"))

    return ds
