import datetime
import math
import warnings

import numpy as np
import xarray as xr
from tqdm import tqdm

from ..core import utils


def ds_rename(ds, waves=False):
    """
    Rename DataArrays within Dataset for EPIC compliance
    """

    varnames = {
        "Pressure": "P_1",
        "pressure": "P_1",
        "Temperature": "Tx_1211",
        "Heading": "Hdg_1215",
        "heading": "Hdg_1215",
        "Pitch": "Ptch_1216",
        "pitch": "Ptch_1216",
        "Roll": "Roll_1217",
        "roll": "Roll_1217",
        "Battery": "Bat_106",
        "batt": "Bat_106",
        "Soundspeed": "SV_80",
    }

    if "Pressure_ac" in ds:
        varnames["Pressure_ac"] = "P_1ac"

    varnames.update(
        {
            "U": "u_1205",
            "V": "v_1206",
            "W": "w_1204",
            "AGC": "AGC_1202",
            "VEL1": "vel1_1277",
            "VEL2": "vel2_1278",
            "VEL3": "vel3_1279",
            "AMP1": "AGC1_1221",
            "AMP2": "AGC2_1222",
            "AMP3": "AGC3_1223",
        }
    )

    for v in varnames:
        if v in ds:
            ds = ds.rename({v: varnames[v]})

    """
    for v in [
        "avgamp1",
        "avgamp2",
        "avgamp3",
        "U",
        "V",
        "W",
        "Depth",
        "water_depth",
        "cellpos",
        #"vel",  # Signature velocity
    ]:
        if v in ds:
            ds = ds.drop_vars(v)
    """

    return ds


def load_cdf(cdf_filename, atmpres=False):
    """
    Load raw .cdf file and, optionally, an atmospheric pressure .cdf file
    """

    ds = xr.load_dataset(cdf_filename)

    # remove units in case we change and we can use larger time steps
    ds.time.encoding.pop("units")

    if atmpres is not False:
        p = xr.load_dataset(atmpres)
        # TODO: check to make sure this data looks OK
        ds["Pressure_ac"] = xr.DataArray(
            ds["Pressure"]
            - p["atmpres"].reindex_like(
                ds["Pressure"], method="nearest", tolerance="5s"
            )
            - p[
                "atmpres"
            ].offset  # need to set a tolerance since we can be off by a couple seconds somewhere
        )

    return ds


def date_parser(year, month, day, hour, minute, second):
    """read csv and parse dates
    https://stackoverflow.com/q/27112591
    """

    return year + "-" + month + "-" + day + " " + hour + ":" + minute + ":" + second


def add_delta_t(ds, waves=False):
    """
    set DELTA_T attribute for EPIC compliance.
    """
    if not waves:
        ds.attrs.update({"DELTA_T": ds.attrs["AQDProfileInterval"]})
    else:
        ds.attrs.update({"DELTA_T": ds.attrs["WaveInterval"]})

    return ds


def make_heading_np(h):
    zero = np.zeros_like(h)
    one = np.ones_like(h)
    return np.array(
        [
            [np.cos(h), np.sin(h), zero],
            [-np.sin(h), np.cos(h), zero],
            [zero, zero, one],
        ]
    )


def make_tilt_np(p, r):
    zero = np.zeros_like(p)
    return np.array(
        [
            [np.cos(p), -np.sin(p) * np.sin(r), -np.cos(r) * np.sin(p)],
            [zero, np.cos(r), -np.sin(r)],
            [np.sin(p), np.sin(r) * np.cos(p), np.cos(p) * np.cos(r)],
        ]
    )


def coord_transform(vel1, vel2, vel3, heading, pitch, roll, T, T_orig, cs, out="ENU"):
    """Perform coordinate transformation"""

    N, M = np.shape(vel1)

    u = np.zeros((N, M))
    v = np.zeros((N, M))
    w = np.zeros((N, M))

    if cs.lower() == out.lower():
        print(
            f"Data already in {cs} coordinates and conversion to {out} requested. Doing nothing."
        )

        u = vel1
        v = vel2
        w = vel3

    elif (
        (cs == "XYZ" and out == "ENU")
        or (cs == "BEAM" and out == "ENU")
        or (cs == "ENU" and out == "BEAM")
    ):
        hh = np.pi * (heading - 90) / 180
        pp = np.pi * pitch / 180
        rr = np.pi * roll / 180

        # compute heading and tilt matrix and reorder dimensions so matrix multiplication works properly
        tilt = make_tilt_np(pp, rr)
        tilt = np.moveaxis(tilt, 2, 0)
        head = make_heading_np(hh)
        head = np.moveaxis(head, 2, 0)

        # compute this inversion so we don't needlessly do it every time in the loop below
        T_orig_inverted = np.linalg.inv(T_orig)
        Rall = head @ tilt @ T

        for i in tqdm(range(N)):
            # select the correct transformation matrix for this burst
            R = Rall[i]

            if cs == "XYZ" and out == "ENU":
                uvw = (
                    R @ T_orig_inverted @ np.array([vel1[i, :], vel2[i, :], vel3[i, :]])
                )
            elif cs == "BEAM" and out == "ENU":
                uvw = R @ np.array([vel1[i, :], vel2[i, :], vel3[i, :]])
            elif cs == "ENU" and out == "BEAM":
                uvw = np.linalg.inv(R) @ np.array([vel1[i, :], vel2[i, :], vel3[i, :]])

            u[i, :] = uvw[0]
            v[i, :] = uvw[1]
            w[i, :] = uvw[2]
    else:
        raise NotImplementedError(
            f"stglib does not currently support input of {cs} and output of {out} coordinates"
        )

    return u, v, w


def swap_bindist_to_depth(ds):
    return ds.swap_dims({"bindist": "depth"})


def set_orientation(VEL, T):
    """
    Create Z variable depending on instrument orientation
    """
    if "Pressure_ac" in VEL:
        presvar = "Pressure_ac"
    else:
        presvar = "Pressure"

    geopotential_datum_name = None

    if "NAVD88_ref" in VEL.attrs or "NAVD88_elevation_ref" in VEL.attrs:
        # if we have NAVD88 elevations of the bed, reference relative to the instrument height in NAVD88
        if "NAVD88_ref" in VEL.attrs:
            elev = VEL.attrs["NAVD88_ref"] + VEL.attrs["transducer_offset_from_bottom"]
        elif "NAVD88_elevation_ref" in VEL.attrs:
            elev = (
                VEL.attrs["NAVD88_elevation_ref"]
                + VEL.attrs["transducer_offset_from_bottom"]
            )
        long_name = "height relative to NAVD88"
        geopotential_datum_name = "NAVD88"
    elif "height_above_geopotential_datum" in VEL.attrs:
        elev = (
            VEL.attrs["height_above_geopotential_datum"]
            + VEL.attrs["transducer_offset_from_bottom"]
        )
        long_name = f"height relative to {VEL.attrs['geopotential_datum_name']}"
        geopotential_datum_name = VEL.attrs["geopotential_datum_name"]
    else:
        # if we don't have NAVD88 elevations, reference to sea-bed elevation
        elev = VEL.attrs["transducer_offset_from_bottom"]
        long_name = "height relative to sea bed"

    T_orig = T.copy()

    if VEL.attrs["orientation"] == "UP":
        print("User instructed that instrument was pointing UP")

        VEL["z"] = xr.DataArray(elev + VEL["bindist"].values, dims="z")
        VEL["depth"] = xr.DataArray(
            np.nanmean(VEL[presvar]) - VEL["bindist"].values, dims="depth"
        )

    if VEL.attrs["orientation"] == "DOWN":
        print("User instructed that instrument was pointing DOWN")
        T[1, :] = -T[1, :]
        T[2, :] = -T[2, :]
        VEL["z"] = xr.DataArray(elev - VEL["bindist"].values, dims="z")
        VEL["depth"] = xr.DataArray(
            np.nanmean(VEL[presvar]) + VEL["bindist"].values, dims="depth"
        )

    VEL["z"].attrs["standard_name"] = "height"
    VEL["z"].attrs["units"] = "m"
    VEL["z"].attrs["positive"] = "up"
    VEL["z"].attrs["axis"] = "Z"
    VEL["z"].attrs["long_name"] = long_name
    if geopotential_datum_name:
        VEL["z"].attrs["geopotential_datum_name"] = geopotential_datum_name

    VEL["depth"].attrs["standard_name"] = "depth"
    VEL["depth"].attrs["units"] = "m"
    VEL["depth"].attrs["positive"] = "down"
    VEL["depth"].attrs["long_name"] = "depth below mean sea level"

    return VEL, T, T_orig


def make_bin_depth(VEL, waves=False):
    """Create bin_depth variable"""

    if "Pressure_ac" in VEL:
        pres = "Pressure_ac"
    else:
        pres = "Pressure"

    if not waves:
        VEL["bin_depth"] = VEL[pres] - VEL["bindist"]
    else:
        VEL["bin_depth"] = VEL[pres].mean(dim="sample") - VEL["bindist"]

    return VEL


def magvar_correct(ds):
    """Correct for magnetic declination at site"""

    if "magnetic_variation_at_site" in ds.attrs:
        magvardeg = ds.attrs["magnetic_variation_at_site"]
    elif "magnetic_variation" in ds.attrs:
        magvardeg = ds.attrs["magnetic_variation"]
    else:
        warnings.warn(
            "No magnetic variation information provided; not correcting compass orientation"
        )
        magvardeg = 0

    if magvardeg != 0:
        histtext = "Rotating heading and horizontal velocities by {} degrees.".format(
            magvardeg
        )

        ds = utils.insert_history(ds, histtext)

        if "Heading" in ds:
            headvar = "Heading"
        elif "Hdg" in ds:
            headvar = "Hdg"
        elif "heading" in ds:
            headvar = "heading"

        ds[headvar].values = ds[headvar].values + magvardeg
        ds[headvar].values[ds[headvar].values >= 360] = (
            ds[headvar].values[ds[headvar].values >= 360] - 360
        )
        ds[headvar].values[ds[headvar].values < 0] = (
            ds[headvar].values[ds[headvar].values < 0] + 360
        )

        if "U" in ds and "V" in ds:
            uvar = "U"
            vvar = "V"
        elif "u_1205" in ds and "v_1206" in ds:
            uvar = "u_1205"
            vvar = "v_1206"

        ds[uvar].values, ds[vvar].values = rotate(
            ds[uvar].values, ds[vvar].values, magvardeg
        )

    return ds


def rotate(u, v, deg):
    rad = np.deg2rad(deg)
    urot = u * np.cos(rad) + v * np.sin(rad)
    vrot = -u * np.sin(rad) + v * np.cos(rad)

    return urot, vrot


def trim_vel(ds, waves=False, data_vars=["U", "V", "W", "AGC"]):
    """Trim velocity data depending on specified method

    Parameters
    ----------
    ds : xarray.Dataset
        The xarray Dataset
    waves: bool, optional
        Flag to determine whether these are waves data. Default False.
    data_vars : array_like
        List of variables to trim. Default ['U', 'V', 'W', 'AGC'].

    Returns
    -------
    xarray.Dataset
        Dataset with trimmed data
    """

    if (
        "trim_method" in ds.attrs
        and ds.attrs["trim_method"].lower() != "none"
        and ds.attrs["trim_method"] is not None
    ):
        if "Pressure_ac" in ds:
            P = ds["Pressure_ac"]
            Ptxt = "atmospherically corrected"
        elif "P_1ac" in ds:
            P = ds["P_1ac"]
            Ptxt = "atmospherically corrected"
        elif "Pressure" in ds:
            # FIXME incorporate press_ ac below
            P = ds["Pressure"]
            Ptxt = "NON-atmospherically corrected"

        if ds.attrs["trim_method"].lower() == "water level":
            for var in data_vars:
                ds[var] = ds[var].where(ds["bindist"] < P)

            histtext = "Trimmed velocity data using {} pressure (water level).".format(
                Ptxt
            )

            ds = utils.insert_history(ds, histtext)

        elif ds.attrs["trim_method"].lower() == "water level sl":
            if "sl_bins" in ds.attrs:
                sl_bins = ds.attrs["sl_bins"]
                sltxt = ds.attrs["sl_bins"]
            elif "sl_bins" not in ds.attrs:
                sl_bins = 0
                sltxt = 0

            for var in data_vars:
                ds[var] = ds[var].where(
                    ds["bindist"]
                    < (P * np.cos(np.deg2rad(ds.attrs["beam_angle"])))
                    - (ds.attrs["bin_size"] * sl_bins)
                )

            histtext = "Trimmed velocity data using {} pressure (water level) and sidelobes (with {} additional surface bins removed).".format(
                Ptxt, sltxt
            )

            ds = utils.insert_history(ds, histtext)

        elif ds.attrs["trim_method"].lower() == "bin range":
            for var in data_vars:
                ds[var] = ds[var].isel(
                    bindist=slice(ds.attrs["good_bins"][0], ds.attrs["good_bins"][1])
                )

            histtext = "Trimmed velocity data using using good_bins of {}.".format(
                ds.attrs["good_bins"]
            )

            ds = utils.insert_history(ds, histtext)

        # find first bin that is all bad values
        # there might be a better way to do this using xarray and named
        # dimensions, but this works for now
        lastbin = np.argmin(np.all(np.isnan(ds[data_vars[0]].values), axis=0) == False)

        if not lastbin == 0:
            # this trims so there are no all-nan rows in the data
            ds = ds.isel(
                bindist=slice(0, lastbin), z=slice(0, lastbin), depth=slice(0, lastbin)
            )

        # TODO: need to add histcomment

        # TODO: add other trim methods
    else:
        print("Did not trim velocity data")

    return ds


def trim_single_bins(ds, var):
    if var + "_trim_single_bins" in ds.attrs:
        print(f"{var}: Trimming single velocity bins")

        singlebins = (~ds[var].isnull()).sum(dim="z") != 1
        ds[var] = ds[var].where(singlebins)

        notetxt = f"Single {var} velocity bins removed. "

        ds = utils.insert_note(ds, var, notetxt)

    return ds


def read_aqd_hdr(basefile):
    """
    Get instrument metadata from .hdr file
    Was formerly readAQDprfHeader.m
    """

    hdrFile = basefile + ".hdr"

    # read in whole file first; shouldn't be too bad since it's small, to
    # check and make sure it's not an HR file.
    # FIXME: this will need to be removed when we start supporting HR
    # check for HR by seeing if extended velocity range is in .hdr
    with open(hdrFile) as f:
        if "Extended velocity range" in f.read():
            raise NotImplementedError(
                "stglib does not currently support Aquadopp HR datasets"
            )

    f = open(hdrFile, "r")
    row = ""

    Instmeta = {}

    while "Hardware configuration" not in row:
        row = f.readline().rstrip()
        if "Profile interval" in row:
            idx = row.find(" sec")
            Instmeta["AQDProfileInterval"] = int(row[38:idx])
        elif "Number of cells" in row:
            Instmeta["AQDNumberOfCells"] = int(row[38:])
        # required here to differentiate from the wave cell size
        elif row.find("Cell size", 0, 9) != -1:
            idx = row.find(" cm")
            Instmeta["AQDCellSize"] = int(row[38:idx])
        elif "Average interval" in row:
            idx = row.find(" sec")
            Instmeta["AQDAverageInterval"] = int(row[38:idx])
        elif "Measurement load" in row:
            idx = row.find(" %")
            Instmeta["AQDMeasurementLoad"] = int(row[38:idx])
        elif "Transmit pulse length" in row:
            idx = row.find(" m")
            Instmeta["AQDTransmitPulseLength"] = float(row[38:idx])
        elif "Blanking distance" in row:
            idx = row.find(" m")
            Instmeta["AQDBlankingDistance"] = float(row[38:idx])
        elif "Compass update rate" in row:
            idx = row.find(" sec")
            Instmeta["AQDCompassUpdateRate"] = int(row[38:idx])
        elif "Wave measurements" in row:
            Instmeta["WaveMeasurements"] = row[38:]
        elif "Wave - Powerlevel" in row:
            Instmeta["WavePower"] = row[38:]
        elif "Wave - Interval" in row:
            idx = row.find(" sec")
            Instmeta["WaveInterval"] = int(row[38:idx])
        elif "Wave - Number of samples" in row:
            Instmeta["WaveNumberOfSamples"] = int(row[38:])
        elif "Wave - Sampling rate" in row:
            Instmeta["WaveSampleRate"] = row[38:]
        elif "Wave - Cell size" in row:
            idx = row.find(" m")
            Instmeta["WaveCellSize"] = float(row[38:idx])
        elif "Analog input 1" in row:
            Instmeta["AQDAnalogInput1"] = row[38:]
        elif "Analog input 2" in row:
            Instmeta["AQDAnalogInput2"] = row[38:]
        elif "Power output" in row:
            Instmeta["AQDAnalogPowerOutput"] = row[38:]
        elif "Powerlevel" in row:
            Instmeta["AQDPowerLevel"] = row[38:]
        elif "Coordinate system" in row:
            Instmeta["AQDCoordinateSystem"] = row[38:]
        elif "Sound speed" in row:
            Instmeta["AQDSoundSpeed"] = row[38:]
        elif "Salinity" in row:
            Instmeta["AQDSalinity"] = row[38:]
        elif "Number of beams" in row:
            Instmeta["AQDNumberOfBeams"] = int(row[38:])
        elif "Number of pings per burst" in row:
            Instmeta["AQDNumberOfPingsPerBurst"] = int(row[38:])
        elif "Software version" in row:
            Instmeta["AQDSoftwareVersion"] = row[38:]
        elif "Deployment name" in row:
            Instmeta["AQDDeploymentName"] = row[38:]
        elif "Deployment time" in row:
            Instmeta["AQDDeploymentTime"] = row[38:]
        elif "Comments" in row:
            Instmeta["AQDComments"] = row[38:]
            # There may be up to three lines of comments, but only if they were added during deployment.
            # These extra lines will be preceded by blanks instead of a field name.
            # After the comments lines are the System lines, which we currenty don't handle,
            # so we can read the next lines and add to Comments if present.
            # Example showing a two-line comment followed by System1 field:
            # Comments                              SP 15916 30 cmab
            #                                       WTS21
            # System1                               19
            for n in range(2):
                row = f.readline().rstrip()
                if len(row) and row[0] == " ":
                    Instmeta["AQDComments"] += "\n"
                    Instmeta["AQDComments"] += row[38:]

    while "Head configuration" not in row:
        row = f.readline().rstrip()
        if "Serial number" in row:
            Instmeta["AQDSerial_Number"] = row[38:]
        elif "Hardware revision" in row:
            Instmeta["AQDHardwareRevision"] = row[38:]
        elif "Revision number" in row:
            Instmeta["AQDRevisionNumber"] = row[38:]
        elif "Recorder size" in row:
            Instmeta["AQDRecorderSize"] = row[38:]
        elif "Firmware version" in row:
            Instmeta["AQDFirmwareVersion"] = row[38:]
        elif "Velocity range" in row:
            Instmeta["AQDVelocityRange"] = row[38:]
        elif "Power output" in row:
            Instmeta["AQDAnalogPowerOutput"] = row[38:]
        elif "Analog input #1 calibration (a0, a1)" in row:
            Instmeta["AQDAnalogInputCal1"] = row[38:]
        elif "Analog input #2 calibration (a0, a1)" in row:
            Instmeta["AQDAnalogInputCal2"] = row[38:]
        elif "Sync signal data out delay" in row:
            Instmeta["AQDSyncOutDelay"] = row[38:]
        elif "Sync signal power down delay" in row:
            Instmeta["AQDSyncPowerDelay"] = row[38:]

    while "Current profile cell center distance from head (m)" not in row:
        row = f.readline().rstrip()
        if "Pressure sensor" in row:
            Instmeta["AQDPressureSensor"] = row[38:]
        elif "Compass" in row:
            Instmeta["AQDCompass"] = row[38:]
        elif "Tilt sensor" in row:
            Instmeta["AQDTilt"] = row[38:]
        elif "Head frequency" in row:
            idx = row.find(" kHz")
            Instmeta["AQDFrequency"] = int(row[38:idx])
        elif "Number of beams" in row:
            Instmeta["AQDNumBeams"] = int(row[38:])
        elif "Serial number" in row:
            Instmeta["AQDHeadSerialNumber"] = row[38:]
        elif "Transformation matrix" in row:
            Instmeta["AQDTransMatrix"] = np.zeros((3, 3))
            for n in np.arange(3):
                Instmeta["AQDTransMatrix"][n, :] = [float(x) for x in row[38:].split()]
                row = f.readline().rstrip()
        elif "Pressure sensor calibration" in row:
            Instmeta["AQDPressureCal"] = row[38:]

    bd = []
    while "Data file format" not in row:
        row = f.readline().rstrip()
        # avoid the header rule line
        if "-" not in row and row != "" and row != "Data file format":
            bd.append(float(row.split()[1]))

    Instmeta["AQDCCD"] = np.array(bd)  # CCD = Cell Center Distance

    # infer some things based on the Aquadopp brochure
    if Instmeta["AQDFrequency"] == 400:
        Instmeta["AQDBeamWidth"] = 3.7
    elif Instmeta["AQDFrequency"] == 600:
        Instmeta["AQDBeamWidth"] = 3.0
    elif Instmeta["AQDFrequency"] == 1000:
        Instmeta["AQDBeamWidth"] = 3.4
    elif Instmeta["AQDFrequency"] == 2000:
        Instmeta["AQDBeamWidth"] = 1.7
    else:
        Instmeta["AQDBeamWidth"] = np.nan

    Instmeta["AQDBeamPattern"] = "convex"
    Instmeta["AQDBeamAngle"] = 25

    f.close()

    return Instmeta


def check_attrs(ds, waves=False, inst_type="AQD"):
    # Add some metadata originally in the run scripts
    ds.attrs["nominal_sensor_depth_note"] = "WATER_DEPTH - " "initial_instrument_height"
    ds.attrs["nominal_sensor_depth"] = (
        ds.attrs["WATER_DEPTH"] - ds.attrs["initial_instrument_height"]
    )
    ds.attrs["transducer_offset_from_bottom"] = ds.attrs["initial_instrument_height"]

    if "initial_instrument_height" not in ds.attrs or np.isnan(
        ds.attrs["initial_instrument_height"]
    ):
        ds.attrs["initial_instrument_height"] = 0

    if inst_type == "AQD":
        ds.attrs["serial_number"] = ds.attrs["AQDSerial_Number"]
        ds.attrs["frequency"] = ds.attrs["AQDFrequency"]
        ds.attrs["beam_width"] = ds.attrs["AQDBeamWidth"]
        ds.attrs["beam_pattern"] = ds.attrs["AQDBeamPattern"]
        ds.attrs["beam_angle"] = ds.attrs["AQDBeamAngle"]
        ds.attrs["salinity_set_by_user"] = ds.attrs["AQDSalinity"]
        ds.attrs["salinity_set_by_user_units"] = "ppt"

        # update metadata from Aquadopp header to CMG standard so that various
        # profilers have the same attribute wording.  Redundant, but necessary
        if not waves:
            ds.attrs["bin_count"] = ds.attrs["AQDNumberOfCells"]
            ds.attrs["bin_size"] = ds.attrs["AQDCellSize"] / 100  # from cm to m
            ds.attrs["blanking_distance"] = ds.attrs["AQDBlankingDistance"]
            # Nortek lists the distance to the center of the first bin as the
            # blanking distance plus one cell size
            ds.attrs["center_first_bin"] = (
                ds.attrs["blanking_distance"] + ds.attrs["bin_size"]
            )  # in m
        else:
            ds.attrs["bin_count"] = 1  # only 1 wave bin
            ds.attrs["bin_size"] = ds.attrs["WaveCellSize"]
            ds.attrs["blanking_distance"] = ds.attrs["AQDBlankingDistance"]
            # TODO: need to set center_first_bin after return in main
            # calling function

        # TODO: Ideally we want to read the error, status, and orientation from
        # the .SEN file, but this requires reading the first good burst, which
        # is unknown at this point.

        ds.attrs["instrument_type"] = "Nortek Aquadopp Profiler"
    elif inst_type == "VEC":
        ds.attrs["serial_number"] = ds.attrs["VECSerialNumber"]
        ds.attrs["frequency"] = ds.attrs["VECFrequency"]
        ds.attrs["instrument_type"] = "Nortek Vector"

    elif inst_type == "SIG":
        ds.attrs["serial_number"] = ds.attrs["SIGSerialNo"]
        ds.attrs["frequency"] = ds.attrs["SIGHeadFrequency"]
        ds.attrs["instrument_type"] = ds.attrs["SIGInstrumentName"]
        if ds.attrs["frequency"] == 1000 or ds.attrs["frequency"] == 500:
            ds.attrs["beam_angle"] = 25
        elif ds.attrs["frequency"] == 250:
            ds.attrs["beam_angle"] = 20
        freq = ds.attrs["frequency"]
        bang = ds.attrs["beam_angle"]
        print(
            f"Signature slant beam acoustic frequency = {freq} with beam_angle = {bang} from vertical"
        )

    return ds


def create_bindist(ds, waves=False):
    """Check instrument orientation and create variables that depend on this"""

    print("Instrument orientation:", ds.attrs["orientation"])
    print(f"{ds.attrs['center_first_bin']=}")
    print(f"{ds.attrs['bin_size']=}")
    print(f"{ds.attrs['bin_count']=}")

    if not waves:
        print(f"Cell center distance = {ds.attrs['AQDCCD']}")
        bindist = ds.attrs["AQDCCD"]
    else:
        bindist = [ds["cellpos"][0]]

    ds["bindist"] = xr.DataArray(bindist, dims=("bindist"), name="bindist")

    return ds


def update_attrs(ds, waves=False):
    """Define dimensions and variables in NetCDF file"""

    ds["latitude"] = xr.DataArray([ds.attrs["latitude"]], dims=("latitude"))
    ds["longitude"] = xr.DataArray([ds.attrs["longitude"]], dims=("longitude"))

    ds["TransMatrix"] = xr.DataArray(ds.attrs["AQDTransMatrix"])
    # Need to remove AQDTransMatrix from attrs for netCDF3 compliance
    ds.attrs.pop("AQDTransMatrix")

    ds["time"].attrs.update(
        {"standard_name": "time", "axis": "T", "long_name": "time (UTC)"}
    )

    ds["latitude"].attrs.update(
        {
            "units": "degree_north",
            "axis": "Y",
            "long_name": "Latitude",
            "standard_name": "latitude",
            "epic_code": 500,
        }
    )

    ds["longitude"].attrs.update(
        {
            "units": "degree_east",
            "axis": "X",
            "long_name": "Longitude",
            "standard_name": "longitude",
            "epic_code": 502,
        }
    )

    # if "position_datum" in ds.attrs:
    #    ds["latitude"].attrs["datum"] = ds.attrs["position_datum"]
    #    ds["longitude"].attrs["datum"] = ds.attrs["position_datum"]

    ds["bindist"].attrs.update(
        {
            "units": "m",
            "long_name": "distance from transducer head",
            "bin_size": ds.attrs["bin_size"],
            "center_first_bin": ds.attrs["center_first_bin"],
            "bin_count": ds.attrs["bin_count"],
            # "transducer_offset_from_bottom": ds.attrs["transducer_offset_from_bottom"],
        }
    )

    ds["Temperature"].attrs.update(
        # {"units": "C", "long_name": "Temperature", "generic_name": "temp"}
        {"units": "C", "long_name": "Temperature"}
    )

    ds["Pressure"].attrs.update(
        {
            "units": "dbar",
            "long_name": "Uncorrected pressure",
            # "generic_name": "press",
            "note": (
                "Raw pressure from instrument, not corrected for changes "
                "in atmospheric pressure"
            ),
        }
    )

    for n in [1, 2, 3]:
        if "VEL" + str(n) in ds:
            ds["VEL" + str(n)].attrs.update(
                {
                    "units": "m s-1",
                    # "transducer_offset_from_bottom": ds.attrs[
                    #    "transducer_offset_from_bottom"
                    # ],
                }
            )
        ds["AMP" + str(n)].attrs.update(
            {
                "long_name": "Beam " + str(n) + " Echo Amplitude",
                "units": "counts",
                # "transducer_offset_from_bottom": ds.attrs[
                #    "transducer_offset_from_bottom"
                # ],
            }
        )

    if not waves:
        veltxt = "current velocity"
    else:
        veltxt = "wave-burst velocity"

    if ds.attrs["AQDCoordinateSystem"] == "ENU":
        ds["U"].attrs["long_name"] = "Eastward " + veltxt
        ds["V"].attrs["long_name"] = "Northward " + veltxt
        ds["W"].attrs["long_name"] = "Vertical " + veltxt
    elif ds.attrs["AQDCoordinateSystem"] == "XYZ":
        ds["X"].attrs["long_name"] = veltxt.capitalize() + " in X Direction"
        ds["Y"].attrs["long_name"] = veltxt.capitalize() + " in Y Direction"
        ds["Z"].attrs["long_name"] = veltxt.capitalize() + " in Z Direction"
    elif ds.attrs["AQDCoordinateSystem"] == "BEAM":
        ds["VEL1"].attrs["long_name"] = "Beam 1 " + veltxt
        ds["VEL2"].attrs["long_name"] = "Beam 2 " + veltxt
        ds["VEL3"].attrs["long_name"] = "Beam 3 " + veltxt

    ds["Battery"].attrs.update({"units": "Volts", "long_name": "Battery Voltage"})

    ds["Pitch"].attrs.update({"units": "degrees", "long_name": "Instrument Pitch"})

    ds["Roll"].attrs.update({"units": "degrees", "long_name": "Instrument Roll"})

    ds["Heading"].attrs.update(
        {
            "units": "degrees",
            "long_name": "Instrument Heading",
            # "datum": "magnetic north",
        }
    )

    # ds["depth"].attrs.update(
    #     {
    #         "units": "m",
    #         "long_name": "mean water depth",
    #         "bin_size": ds.attrs["bin_size"],
    #         "center_first_bin": ds.attrs["center_first_bin"],
    #         "bin_count": ds.attrs["bin_count"],
    #         "transducer_offset_from_bottom": ds.attrs["transducer_offset_from_bottom"],
    #     }
    # )

    ds["TransMatrix"].attrs["long_name"] = "Transformation Matrix " "for this Aquadopp"
    if "burst" in ds:
        ds["burst"].attrs.update({"units": "count", "long_name": "Burst number"})

    return ds


def ds_add_attrs(ds, waves=False, inst_type="AQD"):
    """
    add EPIC and CMG attributes to xarray Dataset
    """

    def add_vel_attributes(vel, dsattrs):
        vel = utils.check_update_attrs(vel, "units", "m s-1")
        if inst_type == "AQD":
            vel.attrs.update(
                {
                    "data_cmnt": (
                        "Velocity in shallowest bin is often suspect and should be used with caution"
                    )
                }
            )

        # TODO: why do we only do trim_method for Water Level SL?
        if (
            "trim_method" in dsattrs
            and dsattrs["trim_method"].lower() == "water level sl"
        ):
            vel.attrs[
                "note"
            ] = "Velocity bins trimmed if out of water or if side lobes intersect sea surface"

    def add_attributes(var, dsattrs):
        if inst_type == "AQD":
            sn = dsattrs["AQDSerial_Number"]
        elif inst_type == "VEC":
            sn = dsattrs["VECSerialNumber"]
        elif inst_type == "SIG":
            sn = dsattrs["SIGSerialNo"]
        """
        var.attrs.update(
            {
                # "serial_number": sn,
                # "initial_instrument_height": dsattrs["initial_instrument_height"],
                # "nominal_instrument_depth": dsattrs["nominal_instrument_depth"],
                # "height_depth_units": "m",
                # "sensor_type": dsattrs["INST_TYPE"],
            }
        )
        """

    # if utils.is_cf(ds):
    #    ds.attrs["featureType"] = "timeSeriesProfile"

    # Update attributes for EPIC and STG compliance
    ds = utils.ds_coord_no_fillvalue(ds)

    ds["time"].attrs.update(
        {"standard_name": "time", "axis": "T", "long_name": "time (UTC)"}
    )

    # ds["depth"].attrs.update(
    #     {
    #         "units": "m",
    #         "long_name": "mean water depth",
    #         "initial_instrument_height": ds.attrs["initial_instrument_height"],
    #         "nominal_instrument_depth": ds.attrs["nominal_instrument_depth"],
    #         "epic_code": 3,
    #     }
    # )

    # if "NAVD88_ref" in ds.attrs:
    #     ds["depth"].attrs["VERT_DATUM"] = "NAVD88"
    #     ds["depth"].attrs["NOTE"] = (
    #         "Computed as platform depth [m NAVD88] "
    #         "- initial_instrument_height - bin "
    #         "distance from transducer"
    #     )

    if "u_1205" in ds:
        ds["u_1205"].attrs.update(
            {
                # "name": "u",
                "long_name": "Eastward Velocity",
                # "generic_name": "u",
                "epic_code": 1205,
            }
        )

    if "v_1206" in ds:
        ds["v_1206"].attrs.update(
            {
                # "name": "v",
                "long_name": "Northward Velocity",
                # "generic_name": "v",
                "epic_code": 1206,
            }
        )

    if "w_1204" in ds:
        ds["w_1204"].attrs.update(
            {
                # "name": "w",
                "long_name": "Vertical Velocity",
                # "generic_name": "w",
                "epic_code": 1204,
            }
        )

    if "w2_1204" in ds:
        ds["w2_1204"].attrs.update(
            {
                # "name": "w",
                "long_name": "Vertical Velocity (2nd)",
                # "generic_name": "w",
                "epic_code": 1204,
            }
        )

    if "AGC_1202" in ds:
        ds["AGC_1202"].attrs.update(
            {
                "units": "counts",
                # "name": "AGC",
                "long_name": "Average Echo Intensity",
                "generic_name": "AGC",
                # "epic_code": 1202,
            }
        )

    if "sample" in ds:
        ds["sample"].encoding["dtype"] = "i4"
        ds["sample"].attrs["long_name"] = "sample number"
        ds["sample"].attrs["units"] = "1"

    if waves:
        if "u_1205" in ds:
            ds["u_1205"].attrs["units"] = "m s-1"

        if "v_1206" in ds:
            ds["v_1206"].attrs["units"] = "m s-1"

        if "w_1204" in ds:
            ds["w_1204"].attrs["units"] = "m s-1"

        if "vel1_1277" in ds:
            ds["vel1_1277"].attrs.update(
                {
                    "units": "m s-1",
                    "long_name": "Beam 1 Velocity",
                    # "generic_name": "vel1",
                    "epic_code": 1277,
                }
            )
        if "vel2_1278" in ds:
            ds["vel2_1278"].attrs.update(
                {
                    "units": "m s-1",
                    "long_name": "Beam 2 Velocity",
                    # "generic_name": "vel2",
                    "epic_code": 1278,
                }
            )
        if "vel3_1279" in ds:
            ds["vel3_1279"].attrs.update(
                {
                    "units": "m s-1",
                    "long_name": "Beam 3 Velocity",
                    # "generic_name": "vel3",
                    "epic_code": 1279,
                }
            )

    if "AGC1_1221" in ds:
        ds["AGC1_1221"].attrs.update(
            {
                "units": "counts",
                "long_name": "Echo Intensity (AGC) Beam 1",
                # "generic_name": "AGC1",
                "epic_code": 1221,
            }
        )

    if "AGC2_1222" in ds:
        ds["AGC2_1222"].attrs.update(
            {
                "units": "counts",
                "long_name": "Echo Intensity (AGC) Beam 2",
                # "generic_name": "AGC2",
                "epic_code": 1222,
            }
        )

    if "AGC3_1223" in ds:
        ds["AGC3_1223"].attrs.update(
            {
                "units": "counts",
                "long_name": "Echo Intensity (AGC) Beam 3",
                # "generic_name": "AGC3",
                "epic_code": 1223,
            }
        )

    if "P_1" in ds:
        ds["P_1"].attrs.update(
            {
                "units": "dbar",
                # "name": "P",
                "long_name": "Uncorrected pressure",
                # "generic_name": "depth",
                "epic_code": 1,
            }
        )  # TODO: is this generic name correct?

    if "P_1ac" in ds:
        ds["P_1ac"].attrs.update(
            # {"units": "dbar", "name": "Pac", "long_name": "Corrected pressure"}
            {"units": "dbar", "long_name": "Corrected pressure"}
        )
        if "P_1ac_note" in ds.attrs:
            ds["P_1ac"].attrs.update({"note": ds.attrs["P_1ac_note"]})

        add_attributes(ds["P_1ac"], ds.attrs)

        histtext = "Atmospheric pressure compensated."

        ds = utils.insert_history(ds, histtext)

    if "bin_depth" in ds:
        ds["bin_depth"].attrs.update(
            # {"units": "m", "name": "bin depth", "long_name": "bin depth"}
            {"units": "m", "long_name": "bin depth"}
        )

        if "P_1ac" in ds:
            if waves:
                ds["bin_depth"].attrs[
                    "note"
                ] = "Actual depth time series of wave burst bin depths. Calculated as corrected pressure (P_1ac) - bindist."
            else:
                ds["bin_depth"].attrs[
                    "note"
                ] = "Actual depth time series of velocity bins. Calculated as corrected pressure (P_1ac) - bindist."
        else:
            ds["bin_depth"].attrs.update(
                {
                    "note": (
                        "Actual depth time series of velocity bins. Calculated "
                        "as pressure (P_1) - bindist."
                    )
                }
            )

    ds["Tx_1211"].attrs.update(
        {
            "units": "C",
            # "name": "Tx",
            "long_name": "Instrument Internal Temperature",
            # "generic_name": "temp",
            "epic_code": 1211,
        }
    )

    ds["Hdg_1215"].attrs.update(
        {
            "units": "degrees",
            # "name": "Hdg",
            "long_name": "Instrument Heading",
            # "generic_name": "hdg",
            "epic_code": 1215,
        }
    )

    if "magnetic_variation_at_site" in ds.attrs:
        ds["Hdg_1215"].attrs["note"] = (
            "Heading is degrees true. Converted "
            "from magnetic with magnetic variation"
            " of %f." % ds.attrs["magnetic_variation_at_site"]
        )
    elif "magnetic_variation" in ds.attrs:
        ds["Hdg_1215"].attrs["note"] = (
            "Heading is degrees true. Converted "
            "from magnetic with magnetic variation"
            " of %f." % ds.attrs["magnetic_variation"]
        )

    ds["Ptch_1216"].attrs.update(
        {
            "units": "degrees",
            # "name": "Ptch",
            "long_name": "Instrument Pitch",
            # "generic_name": "ptch",
            "epic_code": 1216,
        }
    )

    ds["Roll_1217"].attrs.update(
        {
            "units": "degrees",
            # "name": "Roll",
            "long_name": "Instrument Roll",
            # "generic_name": "roll",
            "epic_code": 1217,
        }
    )

    ds["Bat_106"].attrs.update(
        {"units": "V", "long_name": "Battery voltage", "epic_code": 106}
    )

    if "bindist" in ds:
        if inst_type == "AQD":
            blanking_distance = ds.attrs["AQDBlankingDistance"]
        elif inst_type == "SIG":
            blanking_distance = ds.attrs["SIGBurst_BlankingDistance"]
        ds["bindist"].attrs.update(
            {
                "units": "m",
                "long_name": "distance from transducer head",
                "blanking_distance": blanking_distance,
                "note": (
                    "distance is along profile from instrument head to center of bin"
                ),
            }
        )

    if not waves:
        for v in ["AGC_1202", "u_1205", "v_1206", "w_1204"]:
            if v in ds:
                add_attributes(ds[v], ds.attrs)
        for v in ["u_1205", "v_1206", "w_1204", "w2_1204"]:
            if v in ds:
                add_vel_attributes(ds[v], ds.attrs)
    elif waves:
        for v in [
            "vel1_1277",
            "vel2_1278",
            "vel3_1279",
            "AGC1_1221",
            "AGC2_1222",
            "AGC3_1223",
        ]:
            if v in ds:
                add_attributes(ds[v], ds.attrs)

    for v in [
        "P_1",
        "Tx_1211",
        "Hdg_1215",
        "Ptch_1216",
        "Roll_1217",
        "Bat_106",
        "bin_depth",
        "bindist",
    ]:
        if v in ds:
            add_attributes(ds[v], ds.attrs)

    return ds


def check_valid_config_metadata(metadata, inst_type="AQD"):
    vars = [
        "initial_instrument_height",
        "basefile",
        "filename",
        "Conventions",
    ]

    if inst_type == "AQD":
        vars.append("orientation")

    for k in vars:
        if k not in metadata:
            raise KeyError(f"{k} must be defined, most likely in config.yaml")

    if "CF" not in metadata["Conventions"]:
        raise ValueError(
            "Conventions other than a version of the CF Metadata Conventions are not supported"
        )


def apply_wave_coord_output(ds, T, T_orig):
    # Transform coordinates from ENU to BEAM if necessary
    if "wave_coord_output" in ds.attrs:
        if ds.attrs["wave_coord_output"] != "ENU":
            raise NotImplementedError(
                "Only wave_coord_output to ENU is supported at this time"
            )

        if (ds.attrs["wave_coord_output"] == "ENU") and (
            ds.attrs["AQDCoordinateSystem"] == "ENU"
        ):
            print(
                "Requested wave coordinate output to ENU but coordinate system is already ENU. Doing nothing."
            )
            return ds

        histtext = "Converting from {} to {} at user request.".format(
            ds.attrs["AQDCoordinateSystem"],
            ds.attrs["wave_coord_output"],
        )
        print(histtext)
        ds = utils.insert_history(ds, histtext)
        u, v, w = coord_transform(
            ds["VEL1"].values,
            ds["VEL2"].values,
            ds["VEL3"].values,
            ds["Heading"].values,
            ds["Pitch"].values,
            ds["Roll"].values,
            T,
            T_orig,
            ds.attrs["AQDCoordinateSystem"],
            out=ds.attrs["wave_coord_output"],
        )
        ds["U"] = xr.DataArray(u, dims=("time", "sample"))
        ds["V"] = xr.DataArray(v, dims=("time", "sample"))
        ds["W"] = xr.DataArray(w, dims=("time", "sample"))

        ds = ds.drop(["VEL1", "VEL2", "VEL3"])

    return ds
