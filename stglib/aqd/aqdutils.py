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
        "Burst": "burst",
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
            "COR1": "cor1_1285",
            "COR2": "cor2_1286",
            "COR3": "cor3_1287",
            "COR": "cor_avg",
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

    if atmpres is not False and atmpres is not None:
        ds = atmos_correct(ds, atmpres)

    return ds


def atmos_correct(ds, atmpres):
    met = xr.load_dataset(atmpres)
    # need to save attrs before the subtraction, otherwise they are lost
    attrs = ds["Pressure"].attrs
    # need to set a tolerance since we can be off by a couple seconds somewhere
    # TODO is this still needed?
    ds["Pressure_ac"] = xr.DataArray(
        ds["Pressure"]
        - met["atmpres"].reindex_like(ds["Pressure"], method="nearest", tolerance="5s")
        - met["atmpres"].attrs["offset"]
    )
    ds["Pressure_ac"].attrs = attrs

    ds.attrs["atmospheric_pressure_correction_file"] = atmpres
    ds.attrs["atmospheric_pressure_correction_offset_applied"] = met["atmpres"].attrs[
        "offset"
    ]

    histtext = f"Atmospherically corrected using time-series from {atmpres} and offset of {met['atmpres'].offset}"

    ds = utils.insert_history(ds, histtext)

    # Also add it as a note
    ds["Pressure_ac"].attrs["note"] = histtext

    return ds


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


def set_orientation(VEL, T=None, inst_type="AQD"):
    """
    Create Z variable depending on instrument orientation
    """
    if "Pressure_ac" in VEL:
        presvar = "Pressure_ac"
    elif "Pressure" in VEL:
        presvar = "Pressure"
    else:
        presvar = None

    geopotential_datum_name = None

    if "NAVD88_ref" in VEL.attrs or "NAVD88_elevation_ref" in VEL.attrs:
        # if we have NAVD88 elevations of the bed, reference relative to the instrument height in NAVD88
        if "NAVD88_ref" in VEL.attrs:
            navd88_ref = VEL.attrs["NAVD88_ref"]
        elif "NAVD88_elevation_ref" in VEL.attrs:
            navd88_ref = VEL.attrs["NAVD88_elevation_ref"]

        elev = navd88_ref + VEL.attrs["transducer_offset_from_bottom"]
        if "AnalogInput1_height" in VEL.attrs:
            elev_ai1 = navd88_ref + VEL.attrs["AnalogInput1_height"]
        if "AnalogInput2_height" in VEL.attrs:
            elev_ai2 = navd88_ref + VEL.attrs["AnalogInput2_height"]

        long_name = "height relative to NAVD88"
        geopotential_datum_name = "NAVD88"
    elif "height_above_geopotential_datum" in VEL.attrs:
        hagd = VEL.attrs["height_above_geopotential_datum"]
        elev = hagd + VEL.attrs["transducer_offset_from_bottom"]
        if "AnalogInput1_height" in VEL.attrs:
            elev_ai1 = hagd + VEL.attrs["AnalogInput1_height"]
        if "AnalogInput2_height" in VEL.attrs:
            elev_ai2 = hagd + VEL.attrs["AnalogInput2_height"]

        long_name = f"height relative to {VEL.attrs['geopotential_datum_name']}"
        geopotential_datum_name = VEL.attrs["geopotential_datum_name"]
    else:
        # if we don't have NAVD88 elevations, reference to sea-bed elevation
        elev = VEL.attrs["transducer_offset_from_bottom"]
        if "AnalogInput1_height" in VEL.attrs:
            elev_ai1 = VEL.attrs["AnalogInput1_height"]
        if "AnalogInput2_height" in VEL.attrs:
            elev_ai2 = VEL.attrs["AnalogInput2_height"]

        long_name = "height relative to sea bed"

    if inst_type == "AQD":
        T_orig = T.copy()

    if VEL.attrs["orientation"].upper() == "UP":
        print("User instructed that instrument was pointing UP")

        VEL["z"] = xr.DataArray(elev + VEL["bindist"].values, dims="z")
        if presvar:
            VEL["depth"] = xr.DataArray(
                np.nanmean(VEL[presvar]) - VEL["bindist"].values, dims="depth"
            )

    if VEL.attrs["orientation"].upper() == "DOWN":
        print("User instructed that instrument was pointing DOWN")
        if inst_type == "AQD":
            T[1, :] = -T[1, :]
            T[2, :] = -T[2, :]
        VEL["z"] = xr.DataArray(elev - VEL["bindist"].values, dims="z")
        if presvar:
            VEL["depth"] = xr.DataArray(
                np.nanmean(VEL[presvar]) + VEL["bindist"].values, dims="depth"
            )

    if inst_type == "AQD":
        if "AnalogInput1_height" in VEL.attrs:
            VEL["zai1"] = xr.DataArray([elev_ai1], dims="zai1")
        if "AnalogInput2_height" in VEL.attrs:
            VEL["zai2"] = xr.DataArray([elev_ai2], dims="zai2")

    for z in ["z", "zai1", "zai2"]:
        if z not in VEL:
            continue
        VEL[z].attrs["standard_name"] = "height"
        VEL[z].attrs["units"] = "m"
        VEL[z].attrs["positive"] = "up"
        VEL[z].attrs["axis"] = "Z"
        VEL[z].attrs["long_name"] = long_name
        if geopotential_datum_name:
            VEL[z].attrs["geopotential_datum_name"] = geopotential_datum_name

    if "depth" in VEL:
        VEL["depth"].attrs["standard_name"] = "depth"
        VEL["depth"].attrs["units"] = "m"
        VEL["depth"].attrs["positive"] = "down"
        VEL["depth"].attrs["long_name"] = "depth below mean sea level of deployment"

    if inst_type == "AQD":
        return VEL, T, T_orig
    elif inst_type == "RDI":
        return VEL


def make_bin_depth(VEL, waves=False):
    """Create bin_depth variable"""

    if "Pressure_ac" in VEL:
        pres = "Pressure_ac"
    elif "Pressure" in VEL:
        pres = "Pressure"
    else:
        pres = None

    if not pres:
        return VEL

    if "orientation" in VEL.attrs:
        if VEL.attrs["orientation"].lower() == "down":
            if not waves:
                VEL["bin_depth"] = VEL[pres] + VEL["bindist"]
            else:
                VEL["bin_depth"] = VEL[pres].mean(dim="sample") + VEL["bindist"]

        elif VEL.attrs["orientation"].lower() == "up":
            if not waves:
                VEL["bin_depth"] = VEL[pres] - VEL["bindist"]
            else:
                VEL["bin_depth"] = VEL[pres].mean(dim="sample") - VEL["bindist"]

    else:  # if not specified assume up-looking
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

        ds[headvar] = ds[headvar] + magvardeg
        ds[headvar] = ds[headvar] % 360

        uvar = None
        vvar = None
        if "U" in ds and "V" in ds:
            uvar = "U"
            vvar = "V"
        elif "u_1205" in ds and "v_1206" in ds:
            uvar = "u_1205"
            vvar = "v_1206"
        elif "E" in ds and "N" in ds:
            uvar = "E"
            vvar = "N"

        if uvar is not None and vvar is not None:
            uvarattrsbak = ds[uvar].attrs
            vvarattrsbak = ds[vvar].attrs

            ds[uvar], ds[vvar] = rotate(ds[uvar], ds[vvar], magvardeg)

            ds[uvar].attrs = uvarattrsbak
            ds[vvar].attrs = vvarattrsbak

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

    if ds.attrs["orientation"].upper() == "UP":

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
            elif "P_1" in ds:
                P = ds["P_1"]
                Ptxt = "NON-atmospherically corrected"

            if ds.attrs["trim_method"].lower() == "water level":
                for var in data_vars:
                    ds[var] = ds[var].where(ds["bindist"] < P)

                histtext = (
                    "Trimmed velocity data using {} pressure (water level).".format(
                        Ptxt
                    )
                )

                ds = utils.insert_history(ds, histtext)

            elif ds.attrs["trim_method"].lower() == "water level sl":
                if "trim_surf_bins" in ds.attrs:
                    surf_bins = ds.attrs["trim_surf_bins"]
                    histtext = "Trimmed velocity data using {} pressure (water level) and sidelobes (with {} additional surface bins removed).".format(
                        Ptxt, surf_bins
                    )
                else:
                    surf_bins = 0
                    histtext = "Trimmed velocity data using {} pressure (water level) and sidelobes.".format(
                        Ptxt
                    )

                for var in data_vars:
                    ds[var] = ds[var].where(
                        ds["bindist"]
                        < (P * np.cos(np.deg2rad(ds.attrs["beam_angle"])))
                        - (ds.attrs["bin_size"] * surf_bins)
                    )

                ds = utils.insert_history(ds, histtext)

            elif ds.attrs["trim_method"].lower() == "bin range":
                for var in data_vars:
                    ds[var] = ds[var].isel(
                        bindist=slice(
                            ds.attrs["good_bins"][0], ds.attrs["good_bins"][1]
                        )
                    )

                histtext = "Trimmed velocity data using using good_bins of {}.".format(
                    ds.attrs["good_bins"]
                )

                ds = utils.insert_history(ds, histtext)

            # find first bin that is all bad values
            lastbin = (
                ds[data_vars[0]].isnull().all(dim="time").argmax(dim="bindist").values
            )

            if not lastbin == 0:
                # this trims so there are no all-nan rows in the data
                isel = {}
                for v in ["bindist", "depth", "z"]:
                    if v in ds:
                        isel[v] = slice(0, lastbin)

                ds = ds.isel(isel)

            # TODO: need to add histcomment

            # TODO: add other trim methods
        else:
            print("Did not trim velocity data")

    elif ds.attrs["orientation"].upper() == "DOWN":

        if (
            "trim_method" in ds.attrs
            and ds.attrs["trim_method"].lower() != "none"
            and ds.attrs["trim_method"] is not None
        ):
            if "brange" in ds:
                R = ds["brange"]
                Rtxt = "distance to boundary from brange data"
            elif "brange_file" in ds.attrs:
                dsR = xr.open_dataset(ds.attrs["brange_file"], chunks={"time": 300000})
                R = dsR["brange"]
                tstep = (dsR["time"][1] - dsR["time"][0]).values / np.timedelta64(
                    1, "s"
                )
                R = R.reindex_like(ds, method="nearest", tolerance=f"{int(tstep)} s")
                Rtxt = (
                    f"distance to boundary from brange file {ds.attrs['brange_file']}"
                )

            if ds.attrs["trim_method"].lower() == "brange":
                for var in data_vars:
                    ds[var] = ds[var].where(ds["bindist"] < R)

                histtext = "Trimmed velocity data using {}.".format(Rtxt)

                ds = utils.insert_history(ds, histtext)

            elif ds.attrs["trim_method"].lower() == "brange sl":
                if "trim_bottom_bins" in ds.attrs:
                    bot_bins = ds.attrs["trim_bottom_bins"]
                    histtext = "Trimmed velocity data using {} and sidelobes (with {} additional bottom bins removed).".format(
                        Rtxt, bot_bins
                    )
                else:
                    bot_bins = 0
                    histtext = "Trimmed velocity data using {} and sidelobes.".format(
                        Rtxt
                    )

                for var in data_vars:
                    ds[var] = ds[var].where(
                        ds["bindist"]
                        < (R * np.cos(np.deg2rad(ds.attrs["beam_angle"])))
                        - (ds.attrs["bin_size"] * bot_bins)
                    )

                ds = utils.insert_history(ds, histtext)

            elif ds.attrs["trim_method"].lower() == "inst_ht":
                inst_ht = ds.attrs["initial_instrument_height"]
                IHtxt = "initial instrument height"
                for var in data_vars:
                    ds[var] = ds[var].where(ds["bindist"] < inst_ht)

                histtext = "Trimmed velocity data using {}.".format(IHtxt)

                ds = utils.insert_history(ds, histtext)

            elif ds.attrs["trim_method"].lower() == "inst_ht sl":
                inst_ht = ds.attrs["initial_instrument_height"]
                IHtxt = "initial instrument height"
                if "trim_bottom_bins" in ds.attrs:
                    bot_bins = ds.attrs["trim_bottom_bins"]
                    histtext = "Trimmed velocity data using {} and sidelobes (with {} additional bottom bins removed).".format(
                        IHtxt, bot_bins
                    )
                else:
                    bot_bins = 0
                    histtext = "Trimmed velocity data using {} and sidelobes.".format(
                        IHtxt
                    )

                for var in data_vars:
                    ds[var] = ds[var].where(
                        ds["bindist"]
                        < (inst_ht * np.cos(np.deg2rad(ds.attrs["beam_angle"])))
                        - (ds.attrs["bin_size"] * bot_bins)
                    )

                ds = utils.insert_history(ds, histtext)

            elif ds.attrs["trim_method"].lower() == "bin range":
                for var in data_vars:
                    ds[var] = ds[var].isel(
                        bindist=slice(
                            ds.attrs["good_bins"][0], ds.attrs["good_bins"][1]
                        )
                    )

                histtext = "Trimmed velocity data using using good_bins of {}.".format(
                    ds.attrs["good_bins"]
                )

                ds = utils.insert_history(ds, histtext)

            # find first bin that is all bad values
            lastbin = (
                ds[data_vars[0]].isnull().all(dim="time").argmax(dim="bindist").values
            )

            if not lastbin == 0:
                # this trims so there are no all-nan rows in the data
                isel = {}
                for v in ["bindist", "depth", "z"]:
                    if v in ds:
                        isel[v] = slice(0, lastbin)

                ds = ds.isel(isel)

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

    # read in whole file first to see if it's HR or not.
    hr = False
    with open(hdrFile) as f:
        if "Extended velocity range" in f.read():
            hr = True

    f = open(hdrFile)
    row = ""

    Instmeta = {}

    while "User setup" not in row:
        row = f.readline().rstrip()
        if "Number of checksum errors" in row:
            Instmeta["AQDNumberOfChecksumErrors"] = int(row[38:])

    while "Hardware configuration" not in row:
        row = f.readline().rstrip()
        if hr:
            if "Measurement/Burst interval" in row:
                idx = row.find(" sec")
                Instmeta["AQDHRMeasurementBurstInterval"] = int(row[38:idx])
            if "Cell size" in row:
                idx = row.find(" mm")
                Instmeta["AQDHRCellSize"] = int(row[38:idx])
            if "Orientation" in row:
                Instmeta["AQDHROrientation"] = row[38:]
            if "Distance to surface" in row:
                idx = row.find(" m")
                Instmeta["AQDHRDistanceToSurface"] = float(row[38:idx])
            if "Extended velocity range" in row:
                Instmeta["AQDHRExtendedVelocityRange"] = row[38:]
            if "Pulse distance (Lag1)" in row:
                idx = row.find(" m")
                Instmeta["AQDHRPulseLag1"] = float(row[38:idx])
            if "Pulse distance (Lag2)" in row:
                idx = row.find(" m")
                Instmeta["AQDHRPulseLag2"] = float(row[38:idx])
            if "Profile range" in row:
                idx = row.find(" m")
                Instmeta["AQDHRProfileRange"] = float(row[38:idx])
            if "Horizontal velocity range" in row:
                idx = row.find(" m/s")
                Instmeta["AQDHRHorizontalVelRange"] = float(row[38:idx])
            if "Vertical velocity range" in row:
                idx = row.find(" m/s")
                Instmeta["AQDHRVerticalVelRange"] = float(row[38:idx])
            if "Number of cells" in row:
                idx = row.find(" m/s")
                Instmeta["AQDHRNumberOfCells"] = int(row[38:])
            if "Blanking distance" in row:
                idx = row.find(" m")
                Instmeta["AQDHRBlankingDistance"] = float(row[38:idx])
            if "Burst sampling" in row:
                Instmeta["AQDHRBurstSampling"] = row[38:idx]
            if "Samples per burst" in row:
                Instmeta["AQDHRSamplesPerBurst"] = int(row[38:])
            if "Sampling rate" in row:
                Instmeta["AQDHRSamplingRate"] = row[38:]
            if "Powerlevel first ping" in row:
                Instmeta["AQDHRPowerlevelPing1"] = row[38:]
            if "Powerlevel ping 2" in row:
                Instmeta["AQDHRPowerlevelPing2"] = row[38:]
        else:
            if "Profile interval" in row:
                idx = row.find(" sec")
                Instmeta["AQDProfileInterval"] = int(row[38:idx])
            elif "Number of cells" in row:
                Instmeta["AQDNumberOfCells"] = int(row[38:])
            # required here to differentiate from the wave cell size
            elif row.find("Cell size", 0, 9) != -1:
                idx = row.find(" cm")
                Instmeta["AQDCellSize"] = int(row[38:idx])
            elif "Transmit pulse length" in row:
                idx = row.find(" m")
                Instmeta["AQDTransmitPulseLength"] = float(row[38:idx])
            elif "Blanking distance" in row:
                idx = row.find(" m")
                Instmeta["AQDBlankingDistance"] = float(row[38:idx])
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
            elif "Powerlevel" in row:
                Instmeta["AQDPowerLevel"] = row[38:]

        # these are shared between non hr and hr
        if hr:
            shim = "HR"
        else:
            shim = ""
        if "Average interval" in row:
            idx = row.find(" sec")
            Instmeta[f"AQD{shim}AverageInterval"] = int(row[38:idx])
        elif "Measurement load" in row:
            idx = row.find(" %")
            Instmeta[f"AQD{shim}MeasurementLoad"] = int(row[38:idx])
        elif "Compass update rate" in row:
            idx = row.find(" sec")
            Instmeta[f"AQD{shim}CompassUpdateRate"] = int(row[38:idx])
        elif "Analog input 1" in row:
            Instmeta[f"AQD{shim}AnalogInput1"] = row[38:]
        elif "Analog input 2" in row:
            Instmeta[f"AQD{shim}AnalogInput2"] = row[38:]
        elif "Power output" in row:
            Instmeta[f"AQD{shim}AnalogPowerOutput"] = row[38:]
        elif "Coordinate system" in row:
            Instmeta[f"AQD{shim}CoordinateSystem"] = row[38:]
        elif "Sound speed" in row:
            Instmeta[f"AQD{shim}SoundSpeed"] = row[38:]
        elif "Salinity" in row:
            Instmeta[f"AQD{shim}Salinity"] = row[38:]
        elif "Number of beams" in row:
            Instmeta[f"AQD{shim}NumberOfBeams"] = int(row[38:])
        elif "Number of pings per burst" in row:
            Instmeta[f"AQD{shim}NumberOfPingsPerBurst"] = int(row[38:])
        elif "Software version" in row:
            Instmeta[f"AQD{shim}SoftwareVersion"] = row[38:]
        elif "Deployment name" in row:
            Instmeta[f"AQD{shim}DeploymentName"] = row[38:]
        elif "Deployment time" in row:
            Instmeta[f"AQD{shim}DeploymentTime"] = row[38:]
        elif "Wrap mode" in row:
            Instmeta[f"AQD{shim}WrapMode"] = row[38:]
        elif "Comments" in row:
            Instmeta[f"AQD{shim}Comments"] = row[38:]
            # There may be up to three lines of comments, but only if they were added during deployment.
            # These extra lines will be preceded by blanks instead of a field name.
            # After the comments lines are the System lines, which we currently don't handle,
            # so we can read the next lines and add to Comments if present.
            # Example showing a two-line comment followed by System1 field:
            # Comments                              SP 15916 30 cmab
            #                                       WTS21
            # System1                               19
            for n in range(2):
                row = f.readline().rstrip()
                if len(row) and row[0] == " ":
                    Instmeta[f"AQD{shim}Comments"] += "\n"
                    Instmeta[f"AQD{shim}Comments"] += row[38:]

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

    while "Current profile cell center distance" not in row:
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
    bdv = []
    while "Data file format" not in row:
        row = f.readline().rstrip()
        # avoid the header rule line
        if (
            "-" not in row
            and row != ""
            and row != "Data file format"
            and "Beam" not in row
            and "Distances" not in row
        ):
            bd.append(float(row.split()[1]))
            if hr:
                bdv.append(float(row.split()[2]))

    if not hr:
        Instmeta["AQDCCD"] = np.array(bd)  # CCD = Cell Center Distance
    else:
        Instmeta["AQDCCD"] = np.array(bdv)  # CCD = Cell Center Distance
        Instmeta["AQDCCDBEAM"] = np.array(bd)  # CCD = Cell Center Distance

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


def check_attrs(ds, waves=False, hr=False, inst_type="AQD"):
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
        if hr:
            ds.attrs["salinity_set_by_user"] = ds.attrs["AQDHRSalinity"]
        else:
            ds.attrs["salinity_set_by_user"] = ds.attrs["AQDSalinity"]
        ds.attrs["salinity_set_by_user_units"] = "ppt"

        # update metadata from Aquadopp header to CMG standard so that various
        # profilers have the same attribute wording.  Redundant, but necessary
        if not waves and not hr:
            ds.attrs["bin_count"] = ds.attrs["AQDNumberOfCells"]
            ds.attrs["bin_size"] = ds.attrs["AQDCellSize"] / 100  # from cm to m
            ds.attrs["blanking_distance"] = ds.attrs["AQDBlankingDistance"]
            # Nortek lists the distance to the center of the first bin as the
            # blanking distance plus one cell size
            ds.attrs["center_first_bin"] = (
                ds.attrs["blanking_distance"] + ds.attrs["bin_size"]
            )  # in m
        elif hr and not waves:
            ds.attrs["bin_count"] = ds.attrs["AQDHRNumberOfCells"]
            ds.attrs["bin_size"] = ds.attrs["AQDHRCellSize"] / 1000  # from m to m
            ds.attrs["blanking_distance"] = ds.attrs["AQDHRBlankingDistance"]
            # Nortek lists the distance to the center of the first bin as the
            # blanking distance plus one cell size
            ds.attrs["center_first_bin"] = (
                ds.attrs["blanking_distance"] + ds.attrs["bin_size"]
            )  # in m
        elif waves:
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

        # find bin_size and sample_rate attributes
        if (
            ds.attrs["data_type"].upper() == "BURST"
            or ds.attrs["data_type"].upper() == "IBURST"
        ):
            ds.attrs["bin_size"] = ds.attrs["SIGBurst_CellSize"]
            ds.attrs["sample_rate"] = ds.attrs["SIGBurst_SamplingRate"]
        elif (
            ds.attrs["data_type"].upper() == "BURSTHR"
            or ds.attrs["data_type"].upper() == "IBURSTHR"
        ):
            ds.attrs["bin_size"] = ds.attrs["SIGBurstHR_CellSize"]
            ds.attrs["sample_rate"] = ds.attrs["SIGBurst_SamplingRate"]
        elif ds.attrs["data_type"].upper() == "ECHOSOUNDER":
            ds.attrs["bin_size"] = ds.attrs["SIGEchoSounder_CellSize"]
            ds.attrs["sample_rate"] = ds.attrs["SIGBurst_SamplingRate"]
        elif ds.attrs["data_type"].upper() == "AVERAGE":
            ds.attrs["bin_size"] = ds.attrs["SIGAverage_CellSize"]
        elif ds.attrs["data_type"].upper() == "ALT_AVERAGE":
            ds.attrs["bin_size"] = ds.attrs["SIGAlt_Average_CellSize"]

        if (
            ds.attrs["data_type"].upper() == "IBURST"
            or ds.attrs["data_type"].upper() == "IBURSTHR"
            or ds.attrs["data_type"] == "ECHOSOUNDER"
        ):
            ds.attrs["beam_angle"] = 0
            beam_type = "beam 5"
        else:
            if ds.attrs["frequency"] == 1000 or ds.attrs["frequency"] == 500:
                ds.attrs["beam_angle"] = 25
            elif ds.attrs["frequency"] == 250:
                ds.attrs["beam_angle"] = 20
            beam_type = "slant beam"

        if "sample_rate" in ds.attrs and "sample_interval" not in ds.attrs:
            ds.attrs["sample_interval"] = 1 / ds.attrs["sample_rate"]

        freq = ds.attrs["frequency"]
        bang = ds.attrs["beam_angle"]
        print(
            f"Signature {beam_type} acoustic frequency = {freq} with beam_angle = {bang} from vertical"
        )

    elif inst_type == "RDI":
        ds.attrs["serial_number"] = ds.attrs["RDISerialNumber"]
        ds.attrs["instrument_type"] = "Teledyne RDI WorkHorse ADCP"

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


def update_attrs(ds, waves=False, hr=False):
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

    ds["bindist"].attrs.update(
        {
            "units": "m",
            "long_name": "distance from transducer head",
            "bin_size": ds.attrs["bin_size"],
            "center_first_bin": ds.attrs["center_first_bin"],
            "bin_count": ds.attrs["bin_count"],
        }
    )

    ds["Temperature"].attrs.update({"units": "C", "long_name": "Temperature"})

    ds["Pressure"].attrs.update(
        {
            "units": "dbar",
            "long_name": "Uncorrected pressure",
            "note": (
                "Raw pressure from instrument, not corrected for changes "
                "in atmospheric pressure"
            ),
        }
    )

    if "Soundspeed" in ds:
        ds["Soundspeed"].attrs["units"] = "m s-1"

    for n in [1, 2, 3]:
        if "VEL" + str(n) in ds:
            ds["VEL" + str(n)].attrs.update(
                {
                    "units": "m s-1",
                }
            )
        ds["AMP" + str(n)].attrs.update(
            {
                "long_name": "Beam " + str(n) + " Echo Amplitude",
                "units": "counts",
            }
        )

    if not waves:
        veltxt = "current velocity"
    else:
        veltxt = "wave-burst velocity"

    if hr:
        csname = "AQDHRCoordinateSystem"
    else:
        csname = "AQDCoordinateSystem"
    if ds.attrs[csname] == "ENU":
        ds["U"].attrs["long_name"] = "Eastward " + veltxt
        ds["V"].attrs["long_name"] = "Northward " + veltxt
        ds["W"].attrs["long_name"] = "Vertical " + veltxt
    elif ds.attrs[csname] == "XYZ":
        ds["X"].attrs["long_name"] = veltxt.capitalize() + " in X Direction"
        ds["Y"].attrs["long_name"] = veltxt.capitalize() + " in Y Direction"
        ds["Z"].attrs["long_name"] = veltxt.capitalize() + " in Z Direction"
    elif ds.attrs[csname] == "BEAM":
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
        }
    )

    ds["TransMatrix"].attrs["long_name"] = "Transformation Matrix for this Aquadopp"
    if "burst" in ds:
        ds["burst"].attrs.update({"units": "count", "long_name": "Burst number"})

    return ds


def ds_add_attrs(ds, waves=False, hr=False, inst_type="AQD"):
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
            if "trim_surf_bins" in ds.attrs:
                vel.attrs["note"] = (
                    "Velocity bins trimmed if out of water or if side lobes intersect sea surface (with {} additional surface bins removed).".format(
                        ds.attrs["trim_surf_bins"]
                    )
                )
            else:
                vel.attrs["note"] = (
                    "Velocity bins trimmed if out of water or if side lobes intersect sea surface."
                )

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
                "long_name": "Eastward Velocity",
                "epic_code": 1205,
            }
        )

    if "v_1206" in ds:
        ds["v_1206"].attrs.update(
            {
                "long_name": "Northward Velocity",
                "epic_code": 1206,
            }
        )

    if "w_1204" in ds:
        ds["w_1204"].attrs.update(
            {
                "long_name": "Vertical Velocity",
                "epic_code": 1204,
            }
        )

    if "w2_1204" in ds:
        ds["w2_1204"].attrs.update(
            {
                "long_name": "Vertical Velocity (2nd)",
                # "epic_code": 1204,
            }
        )

    if "AGC_1202" in ds:
        ds["AGC_1202"].attrs.update(
            {
                "units": "counts",
                "long_name": "Average Echo Intensity",
                # "generic_name": "AGC",
            }
        )

    if "sample" in ds:
        ds["sample"].encoding["dtype"] = "i4"
        ds["sample"].attrs["long_name"] = "sample number"
        ds["sample"].attrs["units"] = "1"

    if "burst" in ds:
        ds["burst"].encoding["dtype"] = "i4"
        ds["burst"].attrs["long_name"] = "burst number"
        ds["burst"].attrs["units"] = "1"

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
                "epic_code": 1277,
                "standard_name": "radial_sea_water_velocity_away_from_instrument",
            }
        )
    if "vel2_1278" in ds:
        ds["vel2_1278"].attrs.update(
            {
                "units": "m s-1",
                "long_name": "Beam 2 Velocity",
                "epic_code": 1278,
                "standard_name": "radial_sea_water_velocity_away_from_instrument",
            }
        )
    if "vel3_1279" in ds:
        ds["vel3_1279"].attrs.update(
            {
                "units": "m s-1",
                "long_name": "Beam 3 Velocity",
                "epic_code": 1279,
                "standard_name": "radial_sea_water_velocity_away_from_instrument",
            }
        )

    if "agc" in ds:
        ds["agc"].attrs.update(
            {
                "units": "counts",
                "long_name": "Echo Intensity (AGC)",
            }
        )

    if "AGC1_1221" in ds:
        ds["AGC1_1221"].attrs.update(
            {
                "units": "counts",
                "long_name": "Echo Intensity (AGC) Beam 1",
                "epic_code": 1221,
            }
        )

    if "AGC2_1222" in ds:
        ds["AGC2_1222"].attrs.update(
            {
                "units": "counts",
                "long_name": "Echo Intensity (AGC) Beam 2",
                "epic_code": 1222,
            }
        )

    if "AGC3_1223" in ds:
        ds["AGC3_1223"].attrs.update(
            {
                "units": "counts",
                "long_name": "Echo Intensity (AGC) Beam 3",
                "epic_code": 1223,
            }
        )

    if "amp" in ds:
        ds["amp"].attrs.update(
            {
                "units": "Counts",
                "long_name": "Signal Strength (amplitude)",
            }
        )

    if "snr" in ds:
        ds["snr"].attrs.update(
            {
                "units": "dB",
                "long_name": "Signal to Noise Ratio",
            }
        )

    if "SNR1" in ds:
        ds["SNR1"].attrs.update(
            {
                "units": "dB",
                "long_name": "Signal to Noise Ratio Beam 1",
            }
        )

    if "SNR2" in ds:
        ds["SNR2"].attrs.update(
            {
                "units": "dB",
                "long_name": "Signal to Noise Ratio Beam 2",
            }
        )

    if "SNR3" in ds:
        ds["SNR3"].attrs.update(
            {
                "units": "dB",
                "long_name": "Signal to Noise Ratio Beam 3",
            }
        )

    if "cor" in ds:
        ds["cor"].attrs.update(
            {
                "units": "percent",
                "long_name": "Correlation",
            }
        )

    if "cor1_1285" in ds:
        ds["cor1_1285"].attrs.update(
            {
                "units": "percent",
                "long_name": "Correlation Beam 1",
            }
        )

    if "cor2_1286" in ds:
        ds["cor2_1286"].attrs.update(
            {
                "units": "percent",
                "long_name": "Correlation Beam 2",
            }
        )

    if "cor3_1287" in ds:
        ds["cor3_1287"].attrs.update(
            {
                "units": "percent",
                "long_name": "Correlation Beam 3",
            }
        )

    if "P_1" in ds:
        ds["P_1"].attrs.update(
            {
                "units": "dbar",
                "long_name": "Uncorrected pressure",
                "epic_code": 1,
            }
        )

    if "P_1ac" in ds:
        ds["P_1ac"].attrs.update({"units": "dbar", "long_name": "Corrected pressure"})
        if "P_1ac_note" in ds.attrs:
            ds["P_1ac"].attrs.update({"note": ds.attrs["P_1ac_note"]})

    if "bin_depth" in ds:
        ds["bin_depth"].attrs.update({"units": "m", "long_name": "bin depth"})

        if "orientation" in ds.attrs:
            if ds.attrs["orientation"].lower() == "down":
                sign = "+"
            elif ds.attrs["orientation"].lower() == "up":
                sign = "-"
        else:  # if not specified assume up-looking
            sign = "-"

        if "P_1ac" in ds:
            if waves:
                ds["bin_depth"].attrs[
                    "note"
                ] = f"Actual depth time series of wave burst bin depths. Calculated as corrected pressure (P_1ac) {sign} bindist."
            else:
                ds["bin_depth"].attrs[
                    "note"
                ] = f"Actual depth time series of velocity bins. Calculated as corrected pressure (P_1ac) {sign} bindist."
        else:
            ds["bin_depth"].attrs.update(
                {
                    "note": f"Actual depth time series of velocity bins. Calculated as pressure (P_1) {sign} bindist."
                }
            )

    if "SV_80" in ds:
        ds["SV_80"].attrs.update(
            {
                "units": "m s-1",
                "epic_code": 80,
            }
        )

    ds["Tx_1211"].attrs.update(
        {
            "units": "C",
            "long_name": "Instrument Internal Temperature",
            "epic_code": 1211,
        }
    )

    ds["Hdg_1215"].attrs.update(
        {
            "units": "degrees",
            "long_name": "Instrument Heading",
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
            "long_name": "Instrument Pitch",
            "epic_code": 1216,
        }
    )

    ds["Roll_1217"].attrs.update(
        {
            "units": "degrees",
            "long_name": "Instrument Roll",
            "epic_code": 1217,
        }
    )

    ds["Bat_106"].attrs.update(
        {"units": "V", "long_name": "Battery voltage", "epic_code": 106}
    )

    if "brange" in ds:
        ds["brange"].attrs.update(
            {
                "units": "mm",
                "long_name": "Distance from probe to boundary",
                "standard_name": "height_above_sea_floor",
                "note": "Calculated from average of start and end values for burst",
            }
        )

    if "vrange" in ds:
        ds["vrange"].attrs.update(
            {
                "units": "mm",
                "long_name": "Distance from sample volume to boundary",
                "standard_name": "height_above_sea_floor",
                "note": "Calculated from average of start and end values for burst",
            }
        )

    # can apply this to all instruments if needed, but keeping just for vec now
    if inst_type == "VEC":
        ds["orientation"].attrs.update(
            {
                "units": "1",
                "long_name": "instrument orientation",
                "note": "0 = UP; 1 = DOWN",
            }
        )

    if "bindist" in ds and "blanking_distance" not in ds["bindist"].attrs:
        if inst_type == "AQD" and not hr:
            blanking_distance = ds.attrs["AQDBlankingDistance"]
        elif inst_type == "AQD" and hr:
            blanking_distance = ds.attrs["AQDHRBlankingDistance"]
        elif inst_type == "SIG":
            blanking_distance = ds.attrs["SIGBurst_BlankingDistance"]
        elif inst_type == "RDI":
            blanking_distance = ds.attrs["blank"]
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
        for v in ["u_1205", "v_1206", "w_1204", "w2_1204"]:
            if v in ds:
                add_vel_attributes(ds[v], ds.attrs)

    return ds


def check_valid_config_metadata(metadata, inst_type="AQD"):
    varis = [
        "initial_instrument_height",
        "basefile",
        "filename",
        "Conventions",
    ]

    if inst_type in ("AQD", "VEC"):
        varis.append("orientation")

    if inst_type == "VEC":
        varis.append("pressure_sensor_height")
        varis.append("velocity_sample_volume_height")

    for k in varis:
        if k not in metadata:
            raise KeyError(f"{k} must be defined, most likely in YAML config file")

    if "CF" not in metadata["Conventions"]:
        raise ValueError(
            "Conventions other than a version of the CF Metadata Conventions are not supported"
        )

    if "orientation" in metadata:
        metadata["orientation"] = metadata["orientation"].upper()


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


def fill_agc(ds):
    """
    Fill velocity data with a min or max AGC threshold.
    Average AGC (AGC_1202) is used to fill transformed eastward, northward, and upward velocities (u_1205, v_1206, w_1204).
    """

    # lsit velocities to fill by cor threshold(s)
    uvw = ["u_1205", "v_1206", "w_1204", "w2_1204", "vel_b5"]

    if "velocity_agc_min" in ds.attrs:
        for var in uvw:
            if var in ds.data_vars:
                ds[var] = ds[var].where(ds["AGC_1202"] > ds.attrs["velocity_agc_min"])
                notetxt = "Velocity data filled using AGC_1202 minimum threshold of {} counts.".format(
                    ds.attrs["velocity_agc_min"]
                )

                ds = utils.insert_note(ds, var, notetxt)

        ds = utils.insert_history(ds, notetxt)

    if "velocity_agc_max" in ds.attrs:
        for var in uvw:
            if var in ds.data_vars:
                ds[var] = ds[var].where(ds["AGC_1202"] < ds.attrs["velocity_agc_max"])
                notetxt = "Velocity data filled using AGC_1202 maximum threshold of {} counts.".format(
                    ds.attrs["velocity_agc_max"]
                )

                ds = utils.insert_note(ds, var, notetxt)

        ds = utils.insert_history(ds, notetxt)

    return ds


def fill_cor(ds):
    """
    Fill velocity data with a min average correlation threshold. (HR only)
    Average correlation is used to fill transformed eastward, northward, and upward velocities (u_1205, v_1206, w_1204).
    Also include option to fill Average Echo Intensity (AGC) data using average correlation threshold
    """

    # list velocities to fill by agc threshold(s)
    uvw = ["u_1205", "v_1206", "w_1204", "w2_1204", "vel_b5"]
    notetxt = None
    cor_var = None
    if "velocity_cor_min" in ds.attrs:
        for var in uvw:
            if var == "vel_b5":
                cor_var = "cor_b5"
            else:
                cor_var = "cor_avg"
            if var in ds.data_vars:
                ds[var] = ds[var].where(ds[cor_var] > ds.attrs["velocity_cor_min"])
                notetxt = "Velocity data filled using average correlation minimum threshold of {} percent.".format(
                    ds.attrs["velocity_cor_min"]
                )

                ds = utils.insert_note(ds, var, notetxt)

        if notetxt is not None:
            ds = utils.insert_history(ds, notetxt)

    notetxt = None
    cor_var = None
    if "agc_cor_min" in ds.attrs:
        vars = ["AGC_1202", "amp_avg", "amp_b5"]
        for var in vars:
            if var == "amp_b5":
                cor_var = "cor_b5"
            else:
                cor_var = "cor_avg"
            if var in ds.data_vars:
                ds[var] = ds[var].where(ds[cor_var] > ds.attrs["agc_cor_min"])
                notetxt = " Echo data filled using average correlation minimum threshold of {} percent.".format(
                    ds.attrs["agc_cor_min"]
                )

                ds = utils.insert_note(ds, var, notetxt)
        if notetxt is not None:
            ds = utils.insert_history(ds, notetxt)

    return ds


def average_burst(ds):
    # make dictionary of int:dtype pairs, will need to convert back after taking mean
    dint = {}
    intlist = []
    for var in ds.data_vars:
        if ds[var].dtype in (int, np.int64):
            dtypestr = ds[var].dtype
            d = {var: dtypestr}
            dint.update(d)
            intlist.append(var)

    ds = ds.mean(
        "sample", skipna=True, keep_attrs=True
    )  # take mean across 'sample' dim

    for ivar in intlist:
        ds[ivar] = ds[ivar].astype(
            dint[ivar]
        )  # need to retype to int bc np/xarray changes int to float when averaging

    if "burst" in ds:
        ds["burst"].encoding["dtype"] = "i4"

    return ds


def ds_swap_dims(ds):
    # swap vert dim to z or dim specified by vert_dim in config yaml file
    # need to preserve z attrs because swap_dims will remove them

    if "vert_dim" in ds.attrs:
        vdim = ds.attrs["vert_dim"]
        attrsbak = ds[vdim].attrs
        for v in ds.data_vars:
            if "bins" in ds[v].coords:
                ds[v] = ds[v].swap_dims({"bins": vdim})
            elif "bindist" in ds[v].coords:
                ds[v] = ds[v].swap_dims({"bindist": vdim})

        ds[vdim].attrs = attrsbak

        # axis attr set for z in utils.create_z so need to del if other than z
        if ds.attrs["vert_dim"] != "z":
            ds[vdim].attrs["axis"] = "Z"
            del ds["z"].attrs["axis"]

    else:  # set vert dim to z if not specified
        attrsbak = ds["z"].attrs
        for v in ds.data_vars:
            if "bins" in ds[v].coords:
                ds[v] = ds[v].swap_dims({"bins": "z"})
            elif "bindist" in ds[v].coords:
                ds[v] = ds[v].swap_dims({"bindist": "z"})

        ds["z"].attrs = attrsbak

    if "AnalogInput1" in ds:
        ds["AnalogInput1"] = ds["AnalogInput1"].expand_dims("zai1", axis=-1)

    if "AnalogInput2" in ds:
        ds["AnalogInput2"] = ds["AnalogInput2"].expand_dims("zai2", axis=-1)

    return ds


def ds_checksum_check(ds):
    warn = False

    if np.any(ds["Checksum"] == 1):
        warn = True

    if "AQDNumberOfChecksumErrors" in ds.attrs:
        if ds.attrs["AQDNumberOfChecksumErrors"] != 0:
            warn = True

    if warn:
        warnings.warn(
            "Non-zero checksum values found in data. This indicates a failed checksum and potentially bad data. Proceed with caution."
        )

    return ds
