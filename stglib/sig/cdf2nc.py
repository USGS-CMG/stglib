import re
import time

import numpy as np
import xarray as xr
from dask.diagnostics import ProgressBar

from ..aqd import aqdutils
from ..core import filter, qaqc, utils

# import os


def cdf_to_nc(cdf_filename, atmpres=False):
    """
    Load a raw .cdf file and generate a processed .nc file
    """
    print(f"Loading {cdf_filename}")
    # print(os.listdir())
    start_time = time.time()

    ds = xr.open_dataset(cdf_filename, chunks={"time": 200000, "bindist": 48})
    # check for user specified chunks
    if "chunks" in ds.attrs:
        chunksizes = dict(zip(ds.attrs["chunks"][::2], ds.attrs["chunks"][1::2]))
        ds.close()
        # make sure values are type int, needed because they are written out to global attribute in -raw.cdf file they are str
        for key in chunksizes:
            if isinstance(chunksizes[key], str):
                chunksizes[key] = int(chunksizes[key])
        print(f"Using user specified chunksizes = {chunksizes}")

        ds = xr.open_dataset(cdf_filename, chunks=chunksizes)

    end_time = time.time()
    print(f"Finished loading {cdf_filename} in {end_time-start_time:.1f} seconds")

    ds = aqdutils.check_attrs(ds, inst_type="SIG")

    # Add atmospheric pressure offset
    if atmpres is not False:
        ds = aqdutils.atmos_correct(ds, atmpres)

    ds = utils.create_nominal_instrument_depth(ds)

    # print(ds)

    # create Z depending on orientation
    ds = utils.create_z(ds)

    # Clip data to in/out water times or via good_ens
    # Need to clip data after coord transform when using dolfyn
    ds = utils.clip_ds(ds)

    if (
        ds.attrs["data_type"] == "Burst"
        or ds.attrs["data_type"] == "BurstHR"
        or ds.attrs["data_type"] == "Average"
        or ds.attrs["data_type"] == "Alt_Average"
    ):
        # Create separate vel variables first
        ds["U"] = ds["VelEast"]
        ds["V"] = ds["VelNorth"]
        ds["W1"] = ds["VelUp1"]
        ds["W2"] = ds["VelUp2"]

        ds = aqdutils.magvar_correct(ds)
        ds = aqdutils.trim_vel(ds, data_vars=["U", "V", "W1", "W2"])

        # make beam shape data for beam data (vel,amp & cor)
        ds["beam"] = xr.DataArray(range(1, ds["NBeams"][0].values + 1), dims="beam")
        if "VelBeam1" in ds:
            ds["vel"] = xr.concat(
                [ds[f"VelBeam{i}"] for i in range(1, ds["NBeams"][0].values + 1)],
                dim="beam",
            )
        if "CorBeam1" in ds:
            ds["cor"] = xr.concat(
                [ds[f"CorBeam{i}"] for i in range(1, ds["NBeams"][0].values + 1)],
                dim="beam",
            )

        if "AmpBeam1" in ds:
            ds["amp"] = xr.concat(
                [ds[f"AmpBeam{i}"] for i in range(1, ds["NBeams"][0].values + 1)],
                dim="beam",
            )

    ds = aqdutils.make_bin_depth(ds)

    ds = ds_make_tmat(ds)

    ds = ds_make_magaccel_vars(ds)

    ds = ds_make_ahrs_vars(ds)

    # Rename DataArrays for EPIC compliance
    ds = fix_encoding(ds)
    ds = aqdutils.ds_rename(ds)  # for common variables
    ds = ds_drop(ds)
    ds = ds_rename_sig(ds)  # for signature vars not in aqds or vecs
    # swap vert dim to z or user specified in vert_dim
    ds = aqdutils.ds_swap_dims(ds)

    ds = utils.ds_add_lat_lon(ds)

    # Add EPIC and CMG attributes
    ds = aqdutils.ds_add_attrs(ds, inst_type="SIG")  # for common adcp vars
    ds = ds_add_attrs_sig(ds)  # for signature vars

    ds = fix_encoding(ds)

    # Add DELTA_T for EPIC compliance
    # ds = utils.add_delta_t(ds)

    # Add start_time and stop_time attrs
    ds = utils.add_start_stop_time(ds)

    # Add history showing file used
    ds = utils.add_history(ds)

    ds = utils.add_standard_names(ds)  # add common standard names

    ds = drop_unused_dims(ds)

    ds = utils.ds_add_lat_lon(ds)

    ds = utils.ds_coord_no_fillvalue(ds)

    ds = drop_attrs(ds)

    # Add min/max values
    ds = utils.add_min_max(ds)

    # qaqc
    for var in ds.data_vars:
        # need to do this or else a "coordinates" attribute with value of "burst" hangs around
        # ds[var].encoding["coordinates"] = None

        # check for any filtering first
        ds = filter.apply_butter_filt(ds, var)
        ds = filter.apply_med_filt(ds, var)

        ds = qaqc.trim_min(ds, var)
        ds = qaqc.trim_max(ds, var)
        ds = qaqc.trim_min_diff(ds, var)
        ds = qaqc.trim_min_diff_pct(ds, var)
        ds = qaqc.trim_max_diff(ds, var)
        ds = qaqc.trim_maxabs_diff_2d(ds, var)
        ds = qaqc.trim_max_diff_pct(ds, var)
        ds = qaqc.trim_med_diff(ds, var)
        ds = qaqc.trim_med_diff_pct(ds, var)
        ds = qaqc.trim_max_blip(ds, var)
        ds = qaqc.trim_max_blip_pct(ds, var)
        ds = qaqc.trim_bad_ens(ds, var)
        ds = qaqc.trim_bad_ens_indiv(ds, var)
        ds = qaqc.trim_fliers(ds, var)
        ds = qaqc.trim_warmup(ds, var)

    # after check for masking vars by other vars
    for var in ds.data_vars:
        ds = qaqc.trim_mask(ds, var)
        ds = qaqc.trim_mask_expr(ds, var)

    # fill with AGC and Cor threshold
    ds = aqdutils.fill_agc(ds)
    ds = aqdutils.fill_cor(ds)

    # write out nc file by data_type
    if "prefix" in ds.attrs:
        nc_filename = ds.attrs["prefix"] + ds.attrs["filename"]
    else:
        nc_filename = ds.attrs["filename"]

    if ds.attrs["data_type"] == "Burst" or ds.attrs["data_type"] == "BurstHR":
        nc_out = nc_filename + "b-cal.nc"
        print("writing Burst (b) data to netCDF nc file")
        delayed_obj = ds.to_netcdf(nc_out, compute=False)
        with ProgressBar():
            delayed_obj.compute()
        print("Done writing netCDF file", nc_out)

    elif ds.attrs["data_type"] == "IBurst" or ds.attrs["data_type"] == "IBurstHR":
        nc_out = nc_filename + "b5-cal.nc"
        print("writing IBurst (b5) data to netCDF nc file")
        delayed_obj = ds.to_netcdf(nc_out, compute=False)
        with ProgressBar():
            delayed_obj.compute()
        print("Done writing netCDF file", nc_out)

    elif ds.attrs["data_type"] == "EchoSounder":
        nc_out = nc_filename + "e1-cal.nc"
        print("writing Echo1 (echo1) data to netCDF nc file")
        delayed_obj = ds.to_netcdf(nc_out, compute=False)
        with ProgressBar():
            delayed_obj.compute()
        print("Done writing netCDF file", nc_out)

    elif ds.attrs["data_type"] == "Average":
        nc_out = nc_filename + "avg-cal.nc"
        print("writing Average (avg) data to netCDF nc file")
        delayed_obj = ds.to_netcdf(nc_out, compute=False)
        with ProgressBar():
            delayed_obj.compute()
        print("Done writing netCDF file", nc_out)

    elif ds.attrs["data_type"] == "Alt_Average":
        nc_out = nc_filename + "alt-cal.nc"
        print("writing Alt_Average (avg) data to netCDF nc file")
        delayed_obj = ds.to_netcdf(nc_out, compute=False)
        with ProgressBar():
            delayed_obj.compute()
        print("Done writing netCDF file", nc_out)

    utils.check_compliance(nc_out, conventions=ds.attrs["Conventions"])

    end_time = time.time()
    print(
        f"Proceesing cdf2nc for Signature data type {ds.attrs['data_type']} in deployment {ds.attrs['filename']} completed"
    )
    print(f"elapsed time = {end_time-start_time:.1f} seconds")

    return ds


def drop_unused_dims(ds):
    """only keep dims that will be in the final files"""
    thedims = ["bindist"]  # keep bindist needed for wave processing
    for v in ds.data_vars:
        for x in ds[v].dims:
            thedims.append(x)

    for x in ds.dims:
        if x not in thedims:
            ds = ds.drop_vars(x)

    return ds


def drop_attrs(ds):
    """Drop some global attrs"""
    att2del = []
    for k in ds.attrs:
        if re.search("_Beam2xyz$", k):
            att2del.append(k)

    for k in att2del:
        del ds.attrs[k]

    # Get rid of some of the instrument header attributes for other data_types
    dtlist = [
        "Burst",
        "BurstHR",
        "IBurst",
        "IBurstHR",
        "EchoSounder",
        "Average",
        "Alt_Average",
        "Alt_Burst",
    ]
    if ds.attrs["data_type"] in dtlist:
        dtlist.remove(ds.attrs["data_type"])

    rm = []  # initialize list of attrs to be removed
    for k in dtlist:
        for j in ds.attrs:
            if re.match(f"^SIG{k}_", j):
                rm.append(j)
    for k in rm:
        del ds.attrs[k]

    return ds


def ds_drop(ds):
    """
    Drop old DataArrays from Dataset that won't make it into the final .nc file
    """

    todrop = [
        "ExtStatus",
        "NBeams",
        "NCells",
        "PressureSensorTemperature",
        "RTCTemperature",
        "MagnetometerTemperature",
        "Magnetometer",
        "Accelerometer",
        "VEL1",
        "VEL2",
        "VEL3",
        "AMP1",
        "AMP2",
        "AMP3",
        "VelEast",
        "VelNorth",
        "VelUp1",
        "VelUp2",
        "VelSpeed",
        "VelDirection",
        "VelX",
        "VelY",
        "VelZ1",
        "VelZ2",
        "AHRSRotationMatrix",
        "AHRSGyroX",
        "AHRSGyroY",
        "AHRSGyroZ",
        "AHRSQuaternionW",
        "AHRSQuaternionX",
        "AHRSQuaternionY",
        "AHRSQuaternionZ",
        "VelBeam1",
        "VelBeam2",
        "VelBeam3",
        "VelBeam4",
        "AmpBeam1",
        "AmpBeam2",
        "AmpBeam3",
        "AmpBeam4",
        "CorBeam1",
        "CorBeam2",
        "CorBeam3",
        "CorBeam4",
        "AltimeterStatus",
    ]

    if ("AnalogInput1" in ds.attrs) and (ds.attrs["AnalogInput1"].lower() == "true"):
        todrop.remove("AnalogInput1")

    if ("AnalogInput2" in ds.attrs) and (ds.attrs["AnalogInput2"].lower() == "true"):
        todrop.remove("AnalogInput2")

    return ds.drop_vars([t for t in todrop if t in ds.variables])


def ds_rename_sig(ds, waves=False):
    """
    Rename DataArrays within Dataset for compliance
    """
    varnames = {
        "EnsembleCount": "sample",
        "AmbiguityVel": "ambig_vel",
        "Status": "status",
        "Error": "error",
        "TransmitEnergy": "xmit_energy",
        "U": "u_1205",
        "V": "v_1206",
        "W1": "w_1204",
        "W2": "w2_1204",
        "VelSpeed": "CS_300",
        "VelDirection": "CD_310",
        "VelBeam1": "vel1_1277",
        "VelBeam2": "vel1_1278",
        "VelBeam3": "vel1_1279",
        "VelBeam4": "vel1_1280",
        "AmpBeam1": "Sv1_1233",
        "AmpBeam2": "Sv2_1234",
        "AmpBeam3": "Sv3_1235",
        "AmpBeam4": "Sv4_1236",
        "CorBeam1": "cor1_1285",
        "CorBeam2": "cor2_1286",
        "CorBeam3": "cor3_1287",
        "CorBeam4": "cor4_1288",
        "AltimeterDistanceLE": "brangeLE",
        "AltimeterQualityLE": "alt_quality",
        "AltimeterDistanceAST": "brangeAST",
        "AltimeterQualityAST": "ast_quality",
        "AltimeterTimeOffsetAST": "ast_offset_time",
        "AltimeterPressure": "ast_pressure",
        "NominalCorrelation": "cor_nominal",
        "VelBeam5": "vel_b5",
        "AmpBeam5": "amp_b5",
        "CorBeam5": "cor_b5",
        "Echo": "echo_amp",
        "Headingstd": "Hdg_std",
        "Pitchstd": "Ptch_std",
        "Rollstd": "Roll_std",
        "Pressurestd": "Pres_std",
    }

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
        "vel",  # Signature velocity
    ]:
        if v in ds:
            ds = ds.drop_vars(v)
    """
    return ds


def fix_encoding(ds):
    """ensure we don't set dtypes uint for CF compliance"""
    if "units" in ds["time"].encoding:
        ds["time"].encoding.pop("units")

    # use time step to select time encoding
    tstep = ds["time"][1] - ds["time"][0]

    if tstep < np.timedelta64(1, "m"):

        histtext = f"make time encoding to dtype double because tstep {tstep} seconds is < 1 minute, round to milliseconds first"
        ds = utils.insert_history(ds, histtext)

        # round time to milliseconds first
        ds["time"] = ds["time"].dt.round("ms")
        ds["time"].encoding["dtype"] = "double"

    else:
        histtext = f"make time encoding int because tstep {tstep} seconds is >= 1 minute, round time to seconds first"
        ds = utils.insert_history(ds, histtext)

        # round time to seconds if time interval >= 1 minute
        ds["time"] = ds["time"].dt.round("s")

        # check time to make sure it fits in int32, assume seconds for time units
        utils.check_time_fits_in_int32(ds, "time")

        ds["time"].encoding["dtype"] = "i4"

    if "beam" in ds.dims:
        ds["beam"].encoding["dtype"] = "i4"

    for var in ds.data_vars:
        if ds[var].dtype == "uint32" or ds[var].dtype == "uint8":
            ds[var].encoding["dtype"] = "int32"
        if var in ["cor", "cor_b5"]:
            ds[var].encoding["dtype"] = "float32"
        if ds[var].dtype == "float64":
            ds[var].encoding["dtype"] = "float32"

    return ds


def ds_add_attrs_sig(ds):
    """
    add attributes to xarray Dataset
    """

    if "earth" in ds:
        ds["earth"].attrs.update(
            {
                # "units": "1",
                "long_name": "Earth Reference Frame",
            }
        )

    if "earth4" in ds:
        ds["earth4"].attrs.update(
            {
                # "units": "1",
                "long_name": "Earth Reference Frame for 4 beam ADCP",
            }
        )

    if "inst" in ds:
        ds["inst"].attrs.update(
            {
                # "units": "1",
                "long_name": "Instrument Reference Frame",
            }
        )

    if "inst4" in ds:
        ds["inst4"].attrs.update(
            {
                # "units": "1",
                "long_name": "Instrument Reference Frame for 4 beam ADCP",
            }
        )

    if "beam" in ds:
        ds["beam"].attrs.update(
            {
                "units": "1",
                "long_name": "Beam Reference Frame",
            }
        )

    if "q" in ds:
        ds["q"].attrs.update(
            {
                "units": "1",
                "long_name": "Quaternion Vector Components",
            }
        )

    if "orientmat" in ds:
        ds["orientmat"].attrs.update(
            {
                "units": "1",
                "long_name": "Rotation Matrix from AHRS",
            }
        )

    if "gyro" in ds:
        ds["gyro"].attrs.update(
            {
                "units": "rad s-1",
                "long_name": "Angular Velocity from AHRS",
            }
        )

    if "quaternions" in ds:
        ds["quaternions"].attrs.update(
            {
                "units": "1",
                "long_name": "Quaternions",
            }
        )

    if "mag" in ds:
        ds["mag"].attrs.update(
            {
                "units": "uT",
                "long_name": "Magnetometer Data",
            }
        )

    if "accel" in ds:
        ds["accel"].attrs.update(
            {
                "units": "m s-2",
                "long_name": "Vector Acceleration of Instrument",
            }
        )

    if "status" in ds:
        ds["status"].attrs.update(
            {
                "units": "1",
                "long_name": "Status Code",
            }
        )

    if "error" in ds:
        ds["error"].attrs.update(
            {
                "units": "1",
                "long_name": "Error Code",
            }
        )

    if "vel" in ds:
        ds["vel"].attrs.update(
            {
                "units": "m s-1",
                "standard_name": "radial_sea_water_velocity_away_from_instrument",
                "long_name": "Beam Velocity",
            }
        )

    if "vel_b5" in ds:
        ds["vel_b5"].attrs.update(
            {
                "units": "m s-1",
                "standard_name": "radial_sea_water_velocity_away_instrument",
                "long_name": "Beam5 Velocity",
            }
        )

    if "amp" in ds:
        ds["amp"].attrs.update(
            {
                "units": "dB",
                "standard_name": "signal_intensity_from_multibeam_acoustic_doppler_velocity_sensor_in_sea_water",
                "long_name": "Acoustic Signal Amplitude",
            }
        )

    if "amp_b5" in ds:
        ds["amp_b5"].attrs.update(
            {
                "units": "dB",
                "standard_name": "signal_intensity_from_multibeam_acoustic_doppler_velocity_sensor_in_sea_water",
                "long_name": "Beam5 Acoustic Signal Amplitude",
            }
        )

    if "echo_amp" in ds:
        ds["echo_amp"].attrs.update(
            {
                "units": "dB",
                "standard_name": "signal_intensity_from_multibeam_acoustic_doppler_velocity_sensor_in_sea_water",
                "long_name": "Echosounder Signal Amplitude",
            }
        )

    if "cor" in ds:
        ds["cor"].attrs.update(
            {
                "units": "percent",
                "standard_name": "beam_consistency_indicator_from_multibeam_acoustic_doppler_velocity_profiler_in_sea_water",
                "long_name": "Acoustic Signal Correlation",
            }
        )

    if "cor_b5" in ds:
        ds["cor_b5"].attrs.update(
            {
                "units": "percent",
                "standard_name": "beam_consistency_indicator_from_multibeam_acoustic_doppler_velocity_profiler_in_sea_water",
                "long_name": "Beam5 Acoustic Signal Correlation",
            }
        )

    if "cor_nominal" in ds:
        ds["cor_nominal"].attrs.update(
            {
                "units": "percent",
                "long_name": "Nominal Correlation",
            }
        )

    if "SV_80" in ds:
        ds["SV_80"].attrs.update(
            {
                "units": "m s-1",
                "long_name": "Speed of Sound",
            }
        )

    if "ambig_vel" in ds:
        ds["ambig_vel"].attrs.update(
            {
                "units": "m s-1",
                "long_name": "Ambiguity velocity",
            }
        )

    if "brange" in ds:
        ds["brange"].attrs.update(
            {
                "units": "m",
                "standard_name": "altimeter_range",
                "long_name": "Altimeter Range",
            }
        )

    if "alt_quality" in ds:
        ds["alt_quality"] = ds["ast_quality"] / 100
        ds["alt_quality"].encoding["dtype"] = "float32"
        ds["alt_quality"].attrs.update(
            {
                "units": "dB",
                "long_name": "Altimeter Quality",
            }
        )

    if "brangeAST" in ds:
        ds["brangeAST"].attrs.update(
            {
                "units": "m",
                "standard_name": "altimeter_range",
                "long_name": "Acoustic Surface Tracking (AST) Range",
            }
        )

    if "ast_quality" in ds:
        ds["ast_quality"] = ds["ast_quality"] / 100
        ds["ast_quality"].encoding["dtype"] = "float32"
        ds["ast_quality"].attrs.update(
            {
                "units": "dB",
                "long_name": "Acoustic Surface Tracking (AST) Quality",
            }
        )

    if "ast_pressure" in ds:
        ds["ast_pressure"].attrs.update(
            {
                "units": "dbar",
                "long_name": "Acoustic Surface Tracking (AST) Pressure",
            }
        )

    if "ast_offset_time" in ds:
        ds["ast_offset_time"].attrs.update(
            {
                "units": "s",
                "long_name": "Acoustic Surface Tracking (AST) Offset Time",
            }
        )

    if "xmit_energy" in ds:
        ds["xmit_energy"].attrs.update(
            {
                "units": "dB",
                "long_name": "Transmit Energy",
            }
        )

    if "vel1_1277" in ds:
        ds["vel1_1277"].attrs.update(
            {
                "units": "m s-1",
                "long_name": "Beam 1 Velocity",
                "epic_code": 1277,
            }
        )

    if "vel2_1278" in ds:
        ds["vel2_1278"].attrs.update(
            {
                "units": "m s-1",
                "long_name": "Beam 2 Velocity",
                "epic_code": 1278,
            }
        )
    if "vel3_1279" in ds:
        ds["vel3_1279"].attrs.update(
            {
                "units": "m s-1",
                "long_name": "Beam 3 Velocity",
                "epic_code": 1279,
            }
        )

    if "vel4_1280" in ds:
        ds["vel4_1280"].attrs.update(
            {
                "units": "m s-1",
                "long_name": "Beam 4 Velocity",
                "epic_code": 1280,
            }
        )

    if "Sv1_1233" in ds:
        ds["Sv1_1233"].attrs.update(
            {
                "units": "dB",
                "long_name": "Backscatter Strength Beam 1",
                "epic_code": 1233,
            }
        )

    if "Sv2_1234" in ds:
        ds["Sv2_1233"].attrs.update(
            {
                "units": "dB",
                "long_name": "Backscatter Strength Beam 2",
                "epic_code": 1234,
            }
        )

    if "Sv3_1235" in ds:
        ds["Sv3_1235"].attrs.update(
            {
                "units": "dB",
                "long_name": "Backscatter Strength Beam 3",
                "epic_code": 1235,
            }
        )

    if "Sv4_1236" in ds:
        ds["Sv1_1236"].attrs.update(
            {
                "units": "dB",
                "long_name": "Backscatter Strength Beam 4",
                "epic_code": 1236,
            }
        )

    if "Sv_1232" in ds:
        ds["Sv_1232"].attrs.update(
            {
                "units": "dB",
                "long_name": "Mean Backscatter Strength",
                "epic_code": 1232,
            }
        )

    if "cor1_1285" in ds:
        ds["cor1_1285"].attrs.update(
            {
                "units": "percent",
                "long_name": "Beam 1 Correlation",
                "epic_code": 1285,
            }
        )

    if "cor2_1286" in ds:
        ds["cor2_1286"].attrs.update(
            {
                "units": "percent",
                "long_name": "Beam 2 Correlation",
                "epic_code": 1286,
            }
        )

    if "cor3_1286" in ds:
        ds["cor2_1286"].attrs.update(
            {
                "units": "percent",
                "long_name": "Beam 2 Correlation",
                "epic_code": 1286,
            }
        )

    if "TransMatrix" in ds:
        ds["TransMatrix"].attrs.update(
            {
                # "units": "1",
                "long_name": "Instrument Transformation Matrix",
            }
        )

    if "Hdg_std" in ds:
        ds["Hdg_std"].attrs.update(
            {
                "units": "degrees",
                "long_name": "Instrument heading standard deviation",
                "standard_name": "platform_orientation",
                "cell_methods": "time: standard_deviation",
            }
        )

    if "Ptch_std" in ds:
        ds["Ptch_std"].attrs.update(
            {
                "units": "degrees",
                "long_name": "Instrument pitch standard deviation",
                "standard_name": "platform_pitch",
                "cell_methods": "time: standard_deviation",
            }
        )

    if "Roll_std" in ds:
        ds["Roll_std"].attrs.update(
            {
                "units": "degrees",
                "long_name": "Instrument roll standard deviation",
                "standard_name": "platform_roll",
                "cell_methods": "time: standard_deviation",
            }
        )

    if "Pres_std" in ds:
        ds["Pres_std"].attrs.update(
            {
                "units": "dbar",
                "long_name": "Uncorrected pressure standard deviation",
                "standard_name": "sea_water_pressure",
                "cell_methods": "time: standard_deviation",
            }
        )

    return ds


def ds_make_magaccel_vars(ds):
    """
    add magnetometer and accelerometer data to xarray Dataset
    """
    # make mag and accel vars
    if "Magnetometer" in ds:
        if "inst" not in ds.dims:
            ds["inst"] = xr.DataArray(["X", "Y", "Z"], dims="inst")
        ds["mag"] = xr.DataArray(
            ds["Magnetometer"].values.reshape(-1, 3).transpose(), dims=["inst", "time"]
        )
        ds["mag"] = ds["mag"] / 10  # convert milligauss to microtesla

    if "Accelerometer" in ds:
        if "inst" not in ds.dims:
            ds["inst"] = xr.DataArray(["X", "Y", "Z"], dims="inst")
        ds["accel"] = xr.DataArray(
            ds["Accelerometer"].values.reshape(-1, 3).transpose(), dims=["inst", "time"]
        )
        ds["accel"] = ds["accel"] * 9.81  # convert gravity (g) to m s-2

    return ds


def ds_make_ahrs_vars(ds):
    """
    add AHRS vars to xarray Dataset (if they exists)
    """
    # if AHRS will have these
    # make rotation/orientation matrix if exists in dataset
    if "AHRSRotationMatrix" in ds:
        if "earth" not in ds.dims:
            ds["earth"] = xr.DataArray(["E", "N", "U"], dims="earth")
        if "inst" not in ds.dims:
            ds["inst"] = xr.DataArray(["X", "Y", "Z"], dims="inst")
        ds["orientmat"] = xr.DataArray(
            ds["AHRSRotationMatrix"]
            .values.reshape(
                -1,
                3,
                3,
            )
            .transpose((2, 1, 0)),
            dims=["earth", "inst", "time"],
        )

    # make gyro var
    if "AHRSGyroX" in ds:
        if "inst" not in ds.dims:
            ds["inst"] = xr.DataArray(["X", "Y", "Z"], dims="inst")

        ds["gyro"] = xr.concat(
            [ds[f"AHRSGyro{i}"] for i in ds["inst"].values], dim="inst"
        )

        ds["gyro"] = ds["gyro"] * np.pi / 180  # convert to rad s-1

    if "AHRSQuaternionW" in ds:
        if "q" not in ds.dims:
            ds["q"] = xr.DataArray(["W", "X", "Y", "Z"], dims="q")

        ds["quaternions"] = xr.concat(
            [ds[f"AHRSQuaternion{i}"] for i in ds["q"].values], dim="q"
        )

    return ds


def ds_make_tmat(ds):
    """
    add instrument Transformation Matrix to xarray Dataset
    """
    # rename to TransMatrix
    varnames = {"Beam2xyz": "TransMatrix"}
    for v in varnames:
        if v in ds:
            ds = ds.rename({v: varnames[v]})

    # swap dims in Tmat
    if "TransMatrix" in ds.data_vars:
        n = ds["TransMatrix"].shape[0]
        if n == 4:
            if "inst4" not in ds.dims:
                ds["inst4"] = xr.DataArray(["X", "Y", "Z1", "Z2"], dims="inst4")

            if "beam" not in ds.dims:
                ds["beam"] = xr.DataArray(
                    range(1, ds["NBeams"][0].values + 1), dims="beam"
                )

            v = "TransMatrix"
            tdims = ds[v].dims
            ds[v] = ds[v].swap_dims({tdims[0]: "inst4"})
            ds[v] = ds[v].swap_dims({tdims[1]: "beam"})

    return ds
