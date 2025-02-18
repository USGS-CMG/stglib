import warnings

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from stglib.aqd import aqdutils

from . import core
from .core import qaqc, utils


def mat_to_cdf(metadata):
    """
    Process SonTek IQ .mat data to raw .cdf file
    """

    basefile = metadata["basefile"]

    ds = read_iq(basefile + ".mat")

    # write out metadata first, then deal exclusively with xarray attrs
    ds = utils.write_metadata(ds, metadata)

    del metadata

    ds = create_iqbindist(ds)

    ds = utils.ensure_cf(ds)

    # Shift time to middle/handle clock error
    ds = utils.shift_time(ds, ds.attrs["flowSampleDuration"] / 2)

    # configure file
    cdf_filename = ds.attrs["filename"] + "-raw.cdf"

    ds.to_netcdf(cdf_filename, unlimited_dims=["time"])

    print(f"Finished writing data to {cdf_filename}")

    return ds


def cdf_to_nc(cdf_filename, atmpres=False):
    """
    Load a raw .cdf file and generate a processed .nc file
    """

    # Load raw .cdf data
    ds = xr.open_dataset(cdf_filename)

    ds = update_prefixes(ds)

    if atmpres is not False:
        ds = aqdutils.atmos_correct(ds, atmpres)

    # Clip data to in/out water times or via good_ens
    ds = utils.clip_ds(ds)

    ds = vel_to_ms(ds)

    ds = create_iqbindepth(ds)

    ds = create_iqz(ds)

    ds = clean_iq(ds)

    ds = trim_iqvel(ds)

    ds = fill_snr(ds)

    ds = fill_vbper(ds)

    ds = rename_vars(ds)

    # assign min/max:
    ds = utils.add_min_max(ds)

    ds = utils.add_start_stop_time(ds)

    ds = utils.add_delta_t(ds)

    # add lat/lon coordinates
    ds = utils.ds_add_lat_lon(ds)

    # QAQC
    ds = qaqc.call_qaqc(ds)

    ds = fill_velmean(ds)

    ds = utils.create_z(ds)  # added 7/31/2023

    ds = utils.add_standard_names(ds)

    ds = ds_add_attrs(ds)

    # ds = utils.no_p_create_depth(ds) #commented out 7/31/23

    dropvars = [
        "SampleNumber",
        "SampleTime",
        "Volume_Total",
        "Volume_Positive",
        "Volume_Negative",
        "Vel",
        "HorizontalSkew",
        "PressOffsetAdjust",
    ]
    for k in dropvars:
        if k in ds:
            ds = ds.drop(k)

    # add lat/lon coordinates to each variable
    for var in ds.variables:
        if (var not in ds.coords) and ("time" not in var):
            # ds = utils.add_lat_lon(ds, var)
            # cast as float32
            ds = utils.set_var_dtype(ds, var)

    dsflow = ds.copy()
    dsprof = ds.copy()

    dsflow = dsflow.drop([k for k in dsflow if "Profile_" in k])
    dsflow = dsflow.drop(
        ["bin_along", "bin_across"]
    )  # do not need bin dims for the flow data
    dsprof = dsprof.drop([k for k in dsprof if "Profile_" not in k])

    newvars = {}
    for k in dsprof:
        newvars[k] = k.replace("Profile_", "")

    dsprof = dsprof.rename(newvars)

    # Write to .nc file
    print("Writing cleaned/trimmed data to .nc file")

    nc_filename = dsflow.attrs["filename"] + "flow-a.nc"
    dsflow.to_netcdf(
        nc_filename, unlimited_dims=["time"], encoding={"time": {"dtype": "i4"}}
    )
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])
    print("Done writing netCDF file", nc_filename)

    nc_filename = dsprof.attrs["filename"] + "prof-a.nc"
    dsprof.to_netcdf(
        nc_filename, unlimited_dims=["time"], encoding={"time": {"dtype": "i4"}}
    )
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])
    print("Done writing netCDF file", nc_filename)


def read_iq(filnam):
    """Read SonTek IQ data which has been exported as a Matlab .mat file from IQ
    software into an xarray Dataset

    Parameters
    ----------
    filnam : string
        The SonTek .mat filename

    Returns
    -------
    xarray.Dataset
        An xarray Dataset of the IQ data
    """

    iqmat = core.utils.loadmat(filnam)
    # offset = iqmat['FlowSubData_PrfHeader_0_BlankingDistance']
    # beamdist_0 = np.linspace(offset, offset + \
    # 100*iqmat['FlowSubData_PrfHeader_0_CellSize'], 100)
    ds = {}

    ds["time"] = xr.DataArray(
        iqmat["FlowData_SampleTime"],
        attrs={
            "standard_name": "time",
            "axis": "T",
            # per email from SonTek
            "units": "microseconds since 2000-01-01 00:00:00",
            "calendar": "proleptic_gregorian",
        },
        dims="time",
    )

    ds["velbeam"] = xr.DataArray(
        [1, 2, 3, 4],
        dims="velbeam",
        attrs={"long_name": "velocity beam number", "units": "1"},
    )
    ds["beam"] = xr.DataArray(
        [1, 2, 3, 4, 5],
        dims="beam",
        attrs={"long_name": "beam number", "units": "1"},
    )
    # ds['beamdist_0'] = xr.DataArray(beamdist_0, dims='beamdist_0')
    # attrs = {}

    # need to do this because sometimes the flowsubdata and profile data is
    # one burst longer
    timelen = len(ds["time"])

    for k in iqmat:
        if "__" not in k and "FlowSubData" not in k:
            # print(k, np.shape(iqmat[k]))
            if len(np.ravel(iqmat[k])) == len(ds["time"]):
                ds[k] = xr.DataArray(np.ravel(iqmat[k]), dims="time")
                if k in iqmat["Data_Units"]:
                    ds[k].attrs["units"] = iqmat["Data_Units"][k].replace("/s", " s-1")
            elif "_2_" in k or "_3_" in k:
                ds[k] = xr.DataArray(
                    iqmat[k][0:timelen, :], dims=("time", "bin_across")
                )
                if k in iqmat["Data_Units"]:
                    ds[k].attrs["units"] = iqmat["Data_Units"][k].replace("/s", " s-1")
            elif "_0_" in k or "_1_" in k:
                ds[k] = xr.DataArray(iqmat[k][0:timelen, :], dims=("time", "bin_along"))
                if k in iqmat["Data_Units"]:
                    ds[k].attrs["units"] = iqmat["Data_Units"][k].replace("/s", " s-1")
            elif "FlowData_SNR" in k:
                ds[k] = xr.DataArray(iqmat[k][0:timelen, :], dims=("time", "velbeam"))
                ds[k].attrs["units"] = iqmat["Data_Units"][k]
            elif "FlowData_Vel" in k and "XYZ" not in k:
                ds[k] = xr.DataArray(iqmat[k][0:timelen, :], dims=("time", "velbeam"))
                if k in iqmat["Data_Units"]:
                    ds[k].attrs["units"] = iqmat["Data_Units"][k].replace("/s", " s-1")
            elif "FlowData_VelXYZ" in k:
                ds["Vel_X_Center"] = xr.DataArray(iqmat[k][0:timelen, 0], dims=("time"))
                ds["Vel_Z_Center"] = xr.DataArray(iqmat[k][0:timelen, 1], dims=("time"))
                ds["Vel_X_Left"] = xr.DataArray(iqmat[k][0:timelen, 2], dims=("time"))
                ds["Vel_X_Right"] = xr.DataArray(iqmat[k][0:timelen, 3], dims=("time"))
                if k in iqmat["Data_Units"]:
                    xzvars = [
                        "Vel_X_Center",
                        "Vel_Z_Center",
                        "Vel_X_Left",
                        "Vel_X_Right",
                    ]
                    for var in xzvars:
                        ds[var].attrs["units"] = iqmat["Data_Units"][k].replace(
                            "/s", " s-1"
                        )

            elif "FlowData_NoiseLevel" in k:
                ds[k] = xr.DataArray(iqmat[k][0:timelen, :], dims=("time", "beam"))
                if k in iqmat["Data_Units"]:
                    ds[k].attrs["units"] = iqmat["Data_Units"][k].replace("/s", " s-1")

        elif "FlowSubData" in k:
            if "CellSize" in k or "BlankingDistance" in k:
                ds[k] = xr.DataArray(iqmat[k][0:timelen], dims=("time")) / 1000
                if k in iqmat["Data_Units"]:
                    ds[k].attrs["units"] = iqmat["Data_Units"][k]

    ds["bin_along"] = np.arange(ds["Profile_0_Vel"].shape[1])
    ds["bin_across"] = np.arange(ds["Profile_2_Vel"].shape[1])

    ds = xr.Dataset(ds)

    ds.attrs["mean_velocity_equation_type"] = str(
        iqmat["System_IqSetup"]["flowSetup"]["equationType"]
    )

    meanveleq = {
        "1": "Theoretical",
        "2": "Index",
        "3": "Theoretical (Center Beams Only)",
        "4": "Theoretical (Beam 1 Only)",
        "5": "Theoretical (Beam 2 Only)",
    }

    for v in meanveleq:
        if v in ds.attrs["mean_velocity_equation_type"]:
            ds.attrs["mean_velocity_equation_type"] = meanveleq[v]

    if ds.attrs["mean_velocity_equation_type"] == "Index":
        ds.attrs["equation_velocity_type"] = str(
            iqmat["System_IqSetup"]["flowSetup"]["velocityType"]
        )

    if "equation_velocity_type" in ds.attrs:
        veltype = {
            "0": "VelocityXZ (X-Center)",
            "2": "VelocityXZ (Z-Center)",
            "3": "VelocityXZ (X-Left)",
            "4": "VelocityXZ (X-Right)",
            "5": "Average of all Vx",
            "10": "Beam 1 Velocity",
            "11": "Beam 2 Velocity",
            "12": "Beam 3 Velocity",
            "13": "Beam 4 Velocity",
            "14": "VelocityXZ (X-Center, Beam 1 Only)",
            "15": "VelocityXZ (X-Center, Beam 2 Only)",
        }

        for v in veltype:
            if v in ds.attrs["equation_velocity_type"]:
                ds.attrs["equation_velocity_type"] = veltype[v]

    ds.attrs["survey_point_count"] = int(
        iqmat["System_IqSetup"]["flowSetup"]["surveyPointCount"]
    )
    ds.attrs["survey_point_count"] = int(
        iqmat["System_IqSetup"]["flowSetup"]["surveyPointCount"]
    )
    ds.attrs["IQ_location_Y"] = (
        iqmat["System_IqSetup"]["flowSetup"]["instrument_Y"] / 1000
    )
    ds.attrs["IQ_location_Z"] = (
        iqmat["System_IqSetup"]["flowSetup"]["instrument_Z"] / 1000
    )
    ds.attrs["channel_cross_section_Y"] = (
        iqmat["System_IqSetup"]["flowSetup"]["channel_Y"][
            0 : ds.attrs["survey_point_count"]
        ]
    ) / 1000
    ds.attrs["channel_cross_section_Z"] = (
        iqmat["System_IqSetup"]["flowSetup"]["channel_Z"][
            0 : ds.attrs["survey_point_count"]
        ]
    ) / 1000

    for k in iqmat["System_IqSetup"]["basicSetup"]:
        if "spare" not in k:
            ds.attrs[k] = iqmat["System_IqSetup"]["basicSetup"][k]
    for k in iqmat["System_Id"]:
        ds.attrs[k] = iqmat["System_Id"][k]

    if ds.attrs["InstrumentType"] == "IQ":
        ds.attrs["AlongChannelBeamAngle"] = 25
        ds.attrs["AcrossChannelBeamAngle"] = 60
    else:
        print("Check and update beam angle for Sontek IQ instrument")

    return xr.decode_cf(ds)


def create_iqbindist(ds):
    """
    Generate bin distances from transducer along vertical profile for the along (beams 0,1)
    and across (beams 2,3) bins. Cannot make bindist a dim because it changes by sample due
    to bin size changing based on water depth.
    """
    for bm in range(4):
        if bm < 2:
            bdname = "bin_along"
        else:
            bdname = "bin_across"

        r = range(len(ds[bdname]))
        cells = np.zeros(np.shape(ds[f"Profile_{bm}_Vel"]))
        fsdbd = f"FlowSubData_PrfHeader_{bm}_BlankingDistance"
        fsdcs = f"FlowSubData_PrfHeader_{bm}_CellSize"
        for n in range(len(ds["time"])):
            # blanking distance + 1*bin_size = center of first bin
            # blanking distance + N*bin_size = center of N bin
            # due to 0 index, have to add 1 additional bin size to equation below
            # blanking distance + 1*binsize + (bin_along(or bin_across)*bin_size) = center of each bin
            cells[n, :] = (
                ds[fsdbd][n].values + (r * ds[fsdcs][n].values) + (ds[fsdcs][n].values)
            )
        ds[f"Profile_{bm}_bindist"] = xr.DataArray(cells, dims=("time", bdname))

        ds[f"Profile_{bm}_bindist"].attrs.update(
            {
                "units": "m",
                "long_name": "bin (center) distance from transducer",
                "positive": f"{ds.attrs['orientation']}",
                "note": "distance is along vertical profile of transducer",
            }
        )

    return ds


def update_prefixes(ds):
    newvars = {}
    for k in ds:
        if "FlowData" in k:
            newvars[k] = k.replace("FlowData_", "")

        elif "FlowSubData_PrfHeader" in k:
            newvars[k] = k.replace("FlowSubData_PrfHeader_", "Profile_")

    return ds.rename(newvars)


def vel_to_ms(ds):
    """
    Convert velocity data from mm/s to m/s
    """

    for var in ds:
        if "Vel" in var:
            ds[var] = ds[var] / 1000
            ds[var].attrs["units"] = "m s-1"

    return ds


def create_iqbindepth(ds):
    """
    Generate bin depths reltive to pressure.
    """

    if "P_1ac" in ds:
        pres = "P_1ac"
    else:
        pres = "Pressure"

    for bm in range(4):
        if ds.attrs["orientation"].upper() == "UP":
            ds[f"Profile_{bm}_bindepth"] = ds[pres] - ds[f"Profile_{bm}_bindist"]
        elif ds.attrs["orientation"].upper() == "DOWN":
            ds[f"Profile_{bm}_bindepth"] = ds[pres] + ds[f"Profile_{bm}_bindist"]
        else:
            print("Could not create z for bins, specifiy orientation")
        ds[f"Profile_{bm}_bindepth"].attrs.update(
            {
                "units": "m",
                "long_name": "bin(center) depth relative to sea surface",
                "positive": "down",
                "note": "Distance is along vertical profile of transducer",
            }
        )

    return ds


def create_iqz(ds):
    """
    Generate bin heights relative to geopotential datum
    """
    for bm in range(4):
        if "height_above_geopotential_datum" in ds.attrs:
            if ds.attrs["orientation"].upper() == "DOWN":
                if bm < 2:
                    ds[f"Profile_{bm}_z"] = xr.DataArray(
                        ds.attrs["height_above_geopotential_datum"]
                        + ds.attrs["initial_instrument_height"]
                        - ds[f"Profile_{bm}_bindist"].values,
                        dims=("time", "bin_along"),
                    )
                else:
                    ds[f"Profile_{bm}_z"] = xr.DataArray(
                        ds.attrs["height_above_geopotential_datum"]
                        + ds.attrs["initial_instrument_height"]
                        - ds[f"Profile_{bm}_bindist"].values,
                        dims=("time", "bin_across"),
                    )
            elif ds.attrs["orientation"].upper() == "UP":
                if bm < 2:
                    ds[f"Profile_{bm}_z"] = xr.DataArray(
                        ds.attrs["height_above_geopotential_datum"]
                        + ds.attrs["initial_instrument_height"]
                        + ds[f"Profile_{bm}_bindist"].values,
                        dims=("time", "bin_along"),
                    )
                else:
                    ds[f"Profile_{bm}_z"] = xr.DataArray(
                        ds.attrs["height_above_geopotential_datum"]
                        + ds.attrs["initial_instrument_height"]
                        + ds[f"Profile_{bm}_bindist"].values,
                        dims=("time", "bin_across"),
                    )

            else:
                print("Could not create z for bins, specifiy orientation")

        else:
            print(
                "Could not create z for bins, specify height_above_geopotential_datum"
            )

    return ds


def clean_iq(iq):
    """
    Preliminary data cleaning when SNR < 0
    """

    iq["Vel_Mean"].values[iq["Vel_Mean"] < -214748] = np.nan
    iq["Vel"].values[iq["Vel"] == -214748368] = np.nan
    for bm in range(4):
        pr = f"Profile_{bm}_Vel"
        iq[pr].values[iq[pr] == -214748368] = np.nan
        am = f"Profile_{bm}_Amp"
        iq[am].values[iq[am] == 65535] = np.nan
        st = f"Profile_{bm}_VelStd"
        iq[st].values[iq[st] < 0] = np.nan

    return iq


def trim_iqvel(ds):
    """
    Trim velocity data depending on specified method
    """

    if (
        "trim_method" in ds.attrs
        and ds.attrs["trim_method"].lower() != "none"
        and ds.attrs["trim_method"] is not None
    ):
        if "AdjustedPressure" in ds:
            P = ds["AdjustedPressure"]
            Ptxt = "atmospherically corrected"
        elif "P_1ac" in ds:
            P = ds["P_1ac"]
            Ptxt = "atmospherically corrected"
        elif "Pressure" in ds:
            # FIXME incorporate press_ ac below
            P = ds["Pressure"]
            Ptxt = "NON-atmospherically corrected"

        for bm in range(4):  # beams are different angles
            if bm < 2:
                bmangle = ds.attrs["AlongChannelBeamAngle"]
            else:
                bmangle = ds.attrs["AcrossChannelBeamAngle"]

            if ds.attrs["trim_method"].lower() == "water level":
                ds[f"Profile_{bm}_Vel"] = ds[f"Profile_{bm}_Vel"].where(
                    ds[f"Profile_{bm}_bindist"] < P
                )

                histtext = (
                    "Trimmed velocity data using {} pressure (water level).".format(
                        Ptxt
                    )
                )

            elif ds.attrs["trim_method"].lower() == "water level sl":
                ds[f"Profile_{bm}_Vel"] = ds[f"Profile_{bm}_Vel"].where(
                    ds[f"Profile_{bm}_bindist"] < P * np.cos(np.deg2rad(bmangle))
                )

                histtext = "Trimmed velocity data using {} pressure (water level) and sidelobes.".format(
                    Ptxt
                )

        ds = utils.insert_history(ds, histtext)

    else:
        print("Did not trim velocity data")

    return ds


def fill_snr(ds):
    """
    Fill velocity data with corresponding beam snr value threshold
    """
    if "snr_threshold" in ds.attrs:
        Ptxt = str(ds.attrs["snr_threshold"])

        for var in ds:
            if "Vel" and "Profile" in var:
                for bm in range(4):
                    var = f"Profile_{bm}_Vel"
                    ds[var] = ds[var].where(ds.SNR[:, bm] > ds.attrs["snr_threshold"])

            else:
                ds["Vel"] = ds["Vel"].where(ds.SNR > ds.attrs["snr_threshold"])
                ds["Vel_X_Center"] = ds["Vel_X_Center"].where(
                    (ds.SNR[:, 0] > ds.attrs["snr_threshold"])
                    & (ds.SNR[:, 1] > ds.attrs["snr_threshold"])
                )
                ds["Vel_Z_Center"] = ds["Vel_Z_Center"].where(
                    (ds.SNR[:, 0] > ds.attrs["snr_threshold"])
                    & (ds.SNR[:, 1] > ds.attrs["snr_threshold"])
                )
                ds["Vel_X_Left"] = ds["Vel_X_Left"].where(
                    ds.SNR[:, 2] > ds.attrs["snr_threshold"]
                )
                ds["Vel_X_Right"] = ds["Vel_X_Right"].where(
                    ds.SNR[:, 3] > ds.attrs["snr_threshold"]
                )
                ds["Vel_Mean"] = ds["Vel_Mean"].where(
                    (ds.SNR[:, 0] > ds.attrs["snr_threshold"])
                    & (ds.SNR[:, 1] > ds.attrs["snr_threshold"])
                    & (ds.SNR[:, 2] > ds.attrs["snr_threshold"])
                    & (ds.SNR[:, 3] > ds.attrs["snr_threshold"])
                )

            histtext = "Filled velocity data using snr threshold of {} for corresponding beam(s).".format(
                Ptxt
            )

        for var in ds:
            if "Vel" in var:
                ds = utils.insert_note(ds, var, histtext)

        ds = utils.insert_history(ds, histtext)
    else:
        print("Did not fill velocity data using snr threshold")
    return ds


def fill_vbper(ds):
    """
    Fill stage, area, range, profile velocity, and depth data with corresponding vertical beam percent good threshold
    """

    if "vbper_threshold" in ds.attrs:
        Ptxt = str(ds.attrs["vbper_threshold"])

        histtext = "Filling stage, area, range, and D_3 (depth) data using vertical beam percent good threshold threshold of {}.".format(
            Ptxt
        )

        notetxt = "Filled data using vertical beam percent good threshold threshold of {}.".format(
            Ptxt
        )

        varlist = {"Depth", "Stage", "Area", "Range"}

        for k in varlist:
            ds[k] = ds[k].where(ds.VbPercentGood > ds.attrs["vbper_threshold"])

            ds = utils.insert_note(ds, k, notetxt)

        ds = utils.insert_history(ds, histtext)

    else:
        print(
            "Did not fill stage, area, range, and depth data data using vertical beam percent good threshold"
        )

    return ds


def fill_velmean(ds):
    meanvel = "Vel_Mean"
    velvars = {
        "Vel_X_Center",
        "Vel_Z_Center",
        "Vel_X_Left",
        "Vel_X_Right",
        "vel1_1277",
        "vel2_1278",
        "vel3_1279",
        "vel4_1280",
    }

    print("Filling Vel_Mean data using mask of all velocity variables")

    for k in velvars:
        ds[meanvel] = ds[meanvel].where(~ds[k].isnull())

    notetxt = "Filled Vel_Mean data using mask of all velocity variables."
    ds = utils.insert_note(ds, meanvel, notetxt)

    histtext = "Filled Vel_Mean data using mask of all velocity variables."
    ds = utils.insert_history(ds, histtext)

    return ds


def rename_vars(ds):
    # set up dict of instrument -> EPIC variable names

    ds["vel1_1277"] = ds["Vel"].sel(velbeam=1)
    ds["vel2_1278"] = ds["Vel"].sel(velbeam=2)
    ds["vel3_1279"] = ds["Vel"].sel(velbeam=3)
    ds["vel4_1280"] = ds["Vel"].sel(velbeam=4)

    newvars = {}

    varnames = {
        "Batt": "Bat_106",
        "Temp": "T_28",
        "Pitch": "Ptch_1216",
        "Roll": "Roll_1217",
        "Depth": "D_3",
        "Pressure": "P_1",
        "SoundSpeed": "SV_80",
        "Pressure_ac": "P_1ac",
        "Profile_0_Amp": "Profile_AGC1_1221",
        "Profile_0_Vel": "Profile_vel1_1277",
        "Profile_0_VelStd": "Profile_vel1_1277Std",
        "Profile_0_BlankingDistance": "Profile_blanking_distance1",
        "Profile_0_CellSize": "Profile_bin_size1",
        "Profile_0_bindist": "Profile_bindist1",
        "Profile_0_z": "Profile_z1",
        "Profile_0_bindepth": "Profile_bindepth1",
        "Profile_1_Amp": "Profile_AGC2_1222",
        "Profile_1_Vel": "Profile_vel2_1278",
        "Profile_1_VelStd": "Profile_vel2_1278Std",
        "Profile_1_BlankingDistance": "Profile_blanking_distance2",
        "Profile_1_CellSize": "Profile_bin_size2",
        "Profile_1_bindist": "Profile_bindist2",
        "Profile_1_z": "Profile_z2",
        "Profile_1_bindepth": "Profile_bindepth2",
        "Profile_2_Amp": "Profile_AGC3_1223",
        "Profile_2_Vel": "Profile_vel3_1279",
        "Profile_2_VelStd": "Profile_vel3_1279Std",
        "Profile_2_BlankingDistance": "Profile_blanking_distance3",
        "Profile_2_CellSize": "Profile_bin_size3",
        "Profile_2_bindist": "Profile_bindist3",
        "Profile_2_z": "Profile_z3",
        "Profile_2_bindepth": "Profile_bindepth3",
        "Profile_3_Amp": "Profile_AGC4_1224",
        "Profile_3_Vel": "Profile_vel4_1280",
        "Profile_3_VelStd": "Profile_vel4_1280Std",
        "Profile_3_BlankingDistance": "Profile_blanking_distance4",
        "Profile_3_CellSize": "Profile_bin_size4",
        "Profile_3_bindist": "Profile_bindist4",
        "Profile_3_z": "Profile_z4",
        "Profile_3_bindepth": "Profile_bindepth4",
    }
    # check to make sure they exist before trying to rename
    for k in varnames:
        if k in ds:
            newvars[k] = varnames[k]

    return ds.rename(newvars)


def ds_add_attrs(ds):
    attrsnams = ["InstrumentSubType", "InstrumentFriendlyName"]

    for k in attrsnams:
        del ds.attrs[k]

    ds.attrs["serial_number"] = ds.attrs.pop("SerialNumber")
    ds.attrs["instrument_type"] = (
        ds.attrs.pop("InstrumentFamily") + "-" + ds.attrs.pop("InstrumentType")
    )

    if "positive_direction" not in ds.attrs:
        ds.attrs["positive_direction"] = "Not specified"
        warnings.warn(
            "Define positive_direction attribute in yaml file (direction of x arrow on IQ at field site)."
        )

    if "flood_direction" not in ds.attrs:
        ds.attrs["flood_direction"] = "Not specified"
        warnings.warn(
            "Define flood_direction attribute in yaml file (direction of flood velocities at field site)."
        )

    # Update attributes for EPIC and STG compliance
    ds = utils.ds_coord_no_fillvalue(ds)

    ds["time"].attrs.update(
        {"standard_name": "time", "axis": "T", "long_name": "time (UTC)"}
    )

    ds["bin_across"].attrs.update(
        {
            "long_name": "bin number for across channel/skew beams (beams 3 & 4)",
            "units": "1",
            "bin_size": "size of bin may vary by sample. See variables bin_size3 and bin_size4",
            "bin_count": f"{(len(ds.bin_across))}",
            "blanking_distance": "blanking distance may vary by sample. See profile variables blanking_distance3 and blanking_distance4",
            "note": "bin number is along profile from corresponding transducer",
        }
    )

    ds["bin_along"].attrs.update(
        {
            "long_name": "bin number for along channel beams (beams 1 & 2)",
            "units": "1",
            "bin_size": "size of bin may vary by sample. See variables bin_size1 and bin_size2",
            "bin_count": f"{(len(ds.bin_along))}",
            "blanking_distance": "blanking distance may vary by sample. See profile variables blanking_distance1 and blanking_distance2",
            "note": "bin number is along profile from corresponding transducer",
        }
    )

    ds["velbeam"].attrs.update(
        {
            "long_name": "velocity beam number",
            "units": "1",
            "note": "does not include vertical beam (vb, beam 5)",
        }
    )

    ds["beam"].attrs.update(
        {
            "long_name": "beam number",
            "units": "1",
            "note": "includes vertical beam (vb, beam 5)",
        }
    )

    ds["D_3"].attrs.update(
        {
            "long_name": "depth below sea surface",
            "standard_name": "depth",
            "positive": f"{ds.depth.attrs['positive']}",
        }
    )

    d3_note = "Calculated using vertical beam if VbPercentGood is greater than 30% and measured using pressure sesnor if VbPercentGood is less than 30%. Relative to the top of the instrument. See Sontek-IQ Series instrument manual for deatils."
    if "note" in ds["D_3"].attrs:
        ds["D_3"].attrs["note"] = ds["D_3"].attrs["note"] + d3_note
    else:
        ds["D_3"].attrs["note"] = d3_note

    # descriptions from Sontek-IQ Series User's Manual available at
    # http://info.xylem.com/sontek-iq-manual.html
    ds["Stage"].attrs.update(
        {
            "long_name": "Sea surface height (NAVD88)",
            "standard_name": "sea_surface_height_above_geopotential_datum",
            "geopotential_datum_name": "NAVD88",
            "positive": "up",  # positive should always be up regardless of instrument orientation
        }
    )
    ds["Area"].attrs["long_name"] = "Cross-sectional area of user-defined channel"
    ds["Flow"].attrs.update(
        {
            "long_name": "Flow rate (using defined channel geometry)",
            "positive_dir": f"{ds.attrs['positive_direction']}",
            "flood_dir": f"{ds.attrs['flood_direction']}",
        }
    )

    ds["Vel_Mean"].attrs.update(
        {
            "long_name": "Mean velocity (depth-integrated)",
            "positive_dir": f"{ds.attrs['positive_direction']}",
            "flood_dir": f"{ds.attrs['flood_direction']}",
            "mean_velocity_equation_type": f"{ds.attrs['mean_velocity_equation_type']}",
            "mean_velocity_equation_note": "Mean velocity calculation method",
        }
    )

    if "equation_velocity_type" in ds.attrs:
        ds["Vel_Mean"].attrs.update(
            {"equation_velocity_type": f"{ds.attrs['equation_velocity_type']}"}
        )
    ds["Volume_Total"].attrs[
        "long_name"
    ] = "Total water volume (based on all measured flow)"
    ds["Volume_Positive"].attrs[
        "long_name"
    ] = "Total volume of water in the positive downstream direction"
    ds["Volume_Negative"].attrs[
        "long_name"
    ] = "Total volume of water in the negative upstream direction"

    ds["vel1_1277"].attrs.update(
        {
            "long_name": "Beam 1 current velocity",
        }
    )

    ds["vel2_1278"].attrs.update(
        {
            "long_name": "Beam 2 current velocity",
        }
    )

    ds["vel3_1279"].attrs.update(
        {
            "long_name": "Beam 3 current velocity",
        }
    )

    ds["vel4_1280"].attrs.update(
        {
            "long_name": "Beam 4 current velocity",
        }
    )

    ds["Vel_X_Center"].attrs.update(
        {
            "long_name": "X velocity in center of channel (from beams 1 & 2)",
            "positive_dir": f"{ds.attrs['positive_direction']}",
        }
    )
    ds["Vel_Z_Center"].attrs.update(
        {
            "long_name": "Z velocity in center of channel (from beams 1 & 2)",
            "positive_dir": f"{ds.attrs['orientation']}",
        }
    )
    ds["Vel_X_Left"].attrs.update(
        {
            "long_name": "X velocity along left bank (from beam 3)",
            "positive_dir": f"{ds.attrs['positive_direction']}",
        }
    )
    ds["Vel_X_Right"].attrs.update(
        {
            "long_name": "X velocity along right bank (beam 4)",
            "positive_dir": f"{ds.attrs['positive_direction']}",
        }
    )

    ds["VelStd"].attrs["long_name"] = "Velocity standard deviation"
    ds["SNR"].attrs["long_name"] = "Signal-to-noise ratio"
    ds["NoiseLevel"].attrs["long_name"] = "Acoustic noise level"
    ds["Range"].attrs.update(
        {
            "long_name": "distance to sea surface",
            "positive": f"{ds.attrs['orientation']}",
        }
    )
    range_note = "measured using vertical acoustic beam (beam 5)"
    if "note" in ds["Range"].attrs:
        ds["Range"].attrs["note"] = ds["Range"].attrs["note"] + range_note
    else:
        ds["Range"].attrs["note"] = range_note

    ds["T_28"].attrs.update(
        {
            "long_name": "Temperature",
            "epic_code": "28",
            "units": "degree_C",
            "standard_name": "sea_water_temperature",
        }
    )
    ds["P_1"].attrs.update(
        {
            "long_name": "Uncorrected pressure",
            "epic_code": "1",
            "standard_name": "sea_water_pressure",
            "units": "dbar",
        }
    )

    if "P_1ac" in ds.variables:
        ds["P_1ac"].attrs.update(
            {
                "long_name": "Corrected pressure",
                "units": "dbar",
            }
        )

        p1ac_note = "Measurement with atmospheric pressure removed (see SonTek-IQ User's Manual for details)"

        if "note" in ds["P_1ac"].attrs:
            ds["P_1ac"].attrs["note"] = ds["P_1ac"].attrs["note"] + p1ac_note
        else:
            ds["P_1ac"].attrs["note"] = p1ac_note

    ds["Bat_106"].attrs.update({"long_name": "Battery voltage", "epic_code": "106"})
    ds["Ptch_1216"].attrs["long_name"] = "Pitch angle in degrees"

    # to be UDUNITS compatible
    if ds["Ptch_1216"].attrs["units"] == "deg":
        ds["Ptch_1216"].attrs["units"] = "degrees"
    if ds["Roll_1217"].attrs["units"] == "deg":
        ds["Roll_1217"].attrs["units"] = "degrees"
    ds["Roll_1217"].attrs.update(
        {
            "long_name": "Instrument Roll",
            "epic_code": "1217",
        }
    )
    ds["Ptch_1216"].attrs.update(
        {
            "long_name": "Instrument Pitch",
            "epic_code": "1216",
        }
    )
    ds["VbPercentGood"].attrs["long_name"] = "Vertical beam percent good"
    ds["HorizontalSkew"].attrs["long_name"] = "Horizontal skew"
    ds["SystemInWater"].attrs.update(
        {
            "long_name": "Percentage of sample during which instrument was submerged",
            "note": "100% means it was submerged for entire sample",
        }
    )

    # Profile Variables
    for n in range(4):
        bm = n + 1
        agccode = 1221 + n
        velcode = 1277 + n

        ds[f"Profile_AGC{bm}_{agccode}"].attrs.update(
            {"units": "counts", "long_name": f"Echo Intensity (AGC) beam {bm}"}
        )
        ds[f"Profile_vel{bm}_{velcode}Std"].attrs[
            "long_name"
        ] = f"beam {bm} velocity profile standard deviation"
        ds[f"Profile_vel{bm}_{velcode}"].attrs.update(
            {"long_name": f"beam {bm} current velocity"}
        )
        ds[f"Profile_blanking_distance{bm}"].attrs.update(
            {"long_name": f"beam {bm} blanking distance", "units": "m"}
        )
        ds[f"Profile_bin_size{bm}"].attrs.update(
            {"long_name": f"beam {bm} bin size", "units": "m"}
        )
        if "height_above_geopotential_datum" in ds.attrs:
            ds[f"Profile_z{bm}"].attrs.update(
                {
                    "standard_name": "height",
                    "long_name": f"beam {bm} bin height relative to {ds.attrs['geopotential_datum_name']}",
                    "units": "m",
                    "positive": f"{ds.attrs['orientation']}",
                    "axis": "Z",
                }
            )

    return ds
