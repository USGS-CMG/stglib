import re
import warnings

import numpy as np
import pandas as pd
import xarray as xr
from dask.diagnostics import ProgressBar
from tqdm import tqdm

from ..aqd import aqdutils
from ..core import attrs, filter, qaqc, transform, utils


def cdf_to_nc(cdf_filename, atmpres=False):
    """
    Load a raw .cdf file and generate a processed .nc file
    """

    # Load raw .cdf data
    ds = xr.load_dataset(cdf_filename)

    if atmpres is not False:
        ds = aqdutils.atmos_correct(ds, atmpres)

    # Clip data to in/out water times or via good_ens
    ds = utils.clip_ds(ds)

    ds = add_vec_globatts(ds)

    ds = check_orientation(ds)

    # create bindist to pass through utils.create_z
    ds = create_bindist(ds)

    ds = utils.create_z(ds)

    ds = create_analog_z(ds)

    ds = transform.coord_transform(ds)

    ds = aqdutils.magvar_correct(ds)

    ds = combine_vars(ds)

    if ds.attrs["sample_mode"] == "BURST":
        ds = dist_to_boundary(ds)

    ds = scale_analoginput(ds)

    # Drop unused variables
    ds = ds_drop(ds)

    # Rename DataArrays for EPIC compliance
    ds = aqdutils.ds_rename(ds)

    # Shape into burst if not continuous mode
    if ds.attrs["sample_mode"] == "BURST":
        ds = reshape(ds)

    ds = qaqc.drop_vars(ds)

    # Add EPIC and CMG attributes
    ds = aqdutils.ds_add_attrs(ds, inst_type="VEC")

    ds = attrs.ds_add_attrs(ds)

    for var in ds.data_vars:

        ds = utils.var_comment(ds, var)

        # need to do this or else a "coordinates" attribute with value of "burst" hangs around
        ds[var].encoding["coordinates"] = None

        # check if any filtering before other qaqc
        ds = filter.apply_butter_filt(ds, var)
        ds = filter.apply_med_filt(ds, var)
        ds = qaqc.trim_med_diff(ds, var)
        ds = qaqc.trim_med_diff_pct(ds, var)

        ds = qaqc.trim_min(ds, var)
        ds = qaqc.trim_max(ds, var)
        ds = qaqc.trim_maxabs_diff(ds, var)
        ds = qaqc.trim_min_diff(ds, var)
        ds = qaqc.trim_min_diff_pct(ds, var)
        ds = qaqc.trim_max_diff(ds, var)
        ds = qaqc.trim_max_diff_pct(ds, var)
        ds = qaqc.trim_max_blip(ds, var)
        ds = qaqc.trim_max_blip_pct(ds, var)
        ds = qaqc.trim_bad_ens(ds, var)
        ds = qaqc.trim_bad_ens_indiv(ds, var)
        ds = qaqc.trim_fliers(ds, var)
        ds = qaqc.trim_warmup(ds, var)

    ds = fill_snr(ds)

    ds = fill_cor(ds)

    # after check for masking vars by other vars
    for var in ds.data_vars:
        ds = qaqc.trim_mask(ds, var)
        ds = qaqc.trim_mask_expr(ds, var)

    # Add start_time and stop_time attrs
    ds = utils.add_start_stop_time(ds)

    # Add history showing file used
    ds = utils.add_history(ds)

    ds = utils.add_standard_names(ds)

    ds = utils.ds_add_lat_lon(ds)

    ds = utils.add_min_max(ds)

    ds = time_encoding(ds)

    ds = var_encoding(ds)

    ds = utils.ds_coord_no_fillvalue(ds)

    if "chunks" in ds.attrs:
        chunksizes = dict(zip(ds.attrs["chunks"][::2], ds.attrs["chunks"][1::2]))
        for key in chunksizes:
            if isinstance(chunksizes[key], str):
                chunksizes[key] = int(chunksizes[key])
        print(f"Using user specified chunksizes = {chunksizes}")

        ds = ds.chunk(chunks=chunksizes)

    else:
        ds = ds.chunk(chunks="auto")

    create_nc(ds)

    # average files
    if ds.attrs["sample_mode"].upper() == "CONTINUOUS":
        if "average_interval" in ds.attrs:
            ds = fill_time_gaps(ds)
            ds = ds_make_average_shape_mi(ds)
        else:
            return ds

    if "average_duration" in ds.attrs:
        ds.attrs["average_samples_per_burst"] = int(
            ds.attrs["average_duration"] * ds.attrs["sample_rate"]
        )
        ds = ds.isel(sample=slice(0, ds.attrs["average_samples_per_burst"]))

    # before taking mean - need to extract variable that will need to be vector averaged
    va_vars = ["Hdg_1215", "Ptch_1216", "Roll_1217"]

    # get vector averages
    dsVA = utils.make_vector_average_vars(ds, data_vars=va_vars, dim="sample")

    # add 360 degrees if the vector averaging makes heading negative
    dsVA["Hdg_1215"] = dsVA["Hdg_1215"] % 360

    ds = ds.mean(dim="sample", keep_attrs=True)

    # replace vector averaged tilt variables
    for var in dsVA.data_vars:
        ds[var] = dsVA[var]

    if ds.attrs["sample_mode"].upper() == "BURST":
        if "average_duration" in ds.attrs:
            histtext = f"Create burst averaged data product using user specified duration for avergage {ds.attrs['average_duration']} seconds"
        else:
            histtext = "Create burst averaged data product"

    elif ds.attrs["sample_mode"].upper() == "CONTINUOUS":
        if "average_duration" in ds.attrs:
            histtext = f"Create averaged data product from continously sampled data using user specified interval {ds.attrs['average_interval']} seconds and duration {ds.attrs['average_duration']} seconds"
        else:
            histtext = f"Create averaged data product from continuous sampled data using user specified interval {ds.attrs['average_interval']} seconds"

    ds = utils.insert_history(ds, histtext)

    ds = aqdutils.ds_add_attrs(ds, inst_type="VEC")

    ds = utils.add_standard_names(ds)

    ds = utils.add_min_max(ds)

    ds = drop_more_vars(ds)

    ds = drop_dims(ds)

    ds = reorder_dims(ds)

    ds = ds_remove_inst_attrs(ds)

    ds = time_encoding(ds)

    ds = utils.add_delta_t(ds)

    ds = utils.ds_coord_no_fillvalue(ds)

    create_average_nc(ds)

    return ds


def add_vec_globatts(ds):

    ds.attrs["serial_number"] = ds.attrs["VECSerialNumber"]
    ds.attrs["frequency"] = ds.attrs["VECFrequency"]
    ds.attrs["instrument_type"] = "Nortek Vector"
    ds.attrs["sample_rate"] = ds.attrs["VECSamplingRate"]

    if ds.attrs["VECBurstInterval"] != "CONTINUOUS":
        ds.attrs["sample_mode"] = "BURST"
        ds.attrs["samples_per_burst"] = ds.attrs["VECSamplesPerBurst"]
        ds.attrs["burst_interval"] = ds.attrs["VECBurstInterval"]
    elif ds.attrs["VECBurstInterval"] == "CONTINUOUS":
        ds.attrs["sample_mode"] = ds.attrs["VECBurstInterval"]

    return ds


def check_orientation(ds):
    """
    Use orientation variable to verify the vector canister was mounted correctly on mooring for velocity transformations
    """

    # User orientation refers to probe orientation
    # Nortek status code orientation refers to z-axis positive direction
    # See Nortek "The Comprehensive Manual - Velocimeters"
    # section 3.1.7 Orientation of Vector probes

    userorient = ds.attrs["orientation"]

    # last bit of statuscode is orientation, which is saved as variable for vel transformations (inspried by dolfyn)
    sc = str(ds["orientation"].isel(time=int(len(ds["time"]) / 2)).values)
    if sc == "0":
        scname = "UP"
    elif sc == "1":
        scname = "DOWN"
    headtype = ds.attrs["VECHeadSerialNumber"][0:3]

    print(
        f"Instrument reported {headtype} case with orientation status code {sc} -> z-axis positive {scname} at middle of deployment"
    )

    if userorient == "UP":
        print("User instructed probe is pointing UP (sample volume above probe)")
    elif userorient == "DOWN":
        print("User instructed probe is pointing DOWN (sample volume below probe)")
    else:
        raise ValueError("Could not determine instrument orientation from user input")

    flag = False
    if headtype == "VEC":
        if sc == "0" and userorient == "UP":
            flag = True
        elif sc == "1" and userorient == "DOWN":
            flag = True
    elif headtype == "VCH":
        if sc == "0" and userorient == "DOWN":
            flag = True
        elif sc == "1" and userorient == "UP":
            flag = True

    if flag is False:
        print(
            "User-provided orientation matches orientation status code at middle of deployment"
        )
    elif flag is True:
        warnings.warn(
            "User-provided orientation does not match orientation status code at middle of deployment"
        )

        warnings.warn(
            "Incorrect vector orientation will cause erroneous velocity transformations. Check deployment orientation details and vector manual to verify vector was oriented correctly"
        )

    return ds


def create_bindist(ds):
    """
    Create bindist coordinate variable for utils.create_z
    Use initial_instrument_height to get distance from bed to transducer
    According to manual (pp. 18-19), beams cross 0.157m from center transducer, which is in center of measurement volume

    Vector manual: https://support.nortekgroup.com/hc/en-us/articles/360029839351-The-Comprehensive-Manual-Velocimeters

    """

    ds["bindist"] = xr.DataArray(
        [0.157],
        dims="bindist",
        attrs={
            "units": "m",
            "long_name": "distance from transducer head",
            "bin_size": float(ds.attrs["VECSamplingVolume"][:-2]) / 1000,
            "bin_count": 1,
        },
    )

    return ds


def create_analog_z(ds):
    """
    Create z variable for analog inputs depending on instrument orientation
    Mostly taken from ../aqd/aqdutils.py
    """

    geopotential_datum_name = None

    if "NAVD88_ref" in ds.attrs or "NAVD88_elevation_ref" in ds.attrs:
        # if we have NAVD88 elevations of the bed, reference relative to the instrument height in NAVD88
        if "NAVD88_ref" in ds.attrs:
            navd88_ref = ds.attrs["NAVD88_ref"]
        elif "NAVD88_elevation_ref" in ds.attrs:
            navd88_ref = ds.attrs["NAVD88_elevation_ref"]

        # elev = VEL.attrs["NAVD88_ref"] + VEL.attrs["transducer_offset_from_bottom"]
        if "AnalogInput1_height" in ds.attrs:
            elev_ai1 = navd88_ref + ds.attrs["AnalogInput1_height"]
        if "AnalogInput2_height" in ds.attrs:
            elev_ai2 = navd88_ref + ds.attrs["AnalogInput2_height"]

        long_name = "height relative to NAVD88"
        geopotential_datum_name = "NAVD88"
    elif "height_above_geopotential_datum" in ds.attrs:
        # elev = (
        #     VEL.attrs["height_above_geopotential_datum"]
        #     + VEL.attrs["transducer_offset_from_bottom"]
        # )
        hagd = ds.attrs["height_above_geopotential_datum"]
        if "AnalogInput1_height" in ds.attrs:
            elev_ai1 = hagd + ds.attrs["AnalogInput1_height"]
        if "AnalogInput2_height" in ds.attrs:
            elev_ai2 = hagd + ds.attrs["AnalogInput2_height"]

        long_name = f"height relative to {ds.attrs['geopotential_datum_name']}"
        geopotential_datum_name = ds.attrs["geopotential_datum_name"]
    else:
        # if we don't have NAVD88 elevations, reference to sea-bed elevation
        # elev = VEL.attrs["transducer_offset_from_bottom"]
        if "AnalogInput1_height" in ds.attrs:
            elev_ai1 = ds.attrs["AnalogInput1_height"]
        if "AnalogInput2_height" in ds.attrs:
            elev_ai2 = ds.attrs["AnalogInput2_height"]

        long_name = "height relative to sea bed"

    if "AnalogInput1_height" in ds.attrs:
        ds["zai1"] = xr.DataArray([elev_ai1], dims="zai1")
    if "AnalogInput2_height" in ds.attrs:
        ds["zai2"] = xr.DataArray([elev_ai2], dims="zai2")

    lnshim = {
        "zai1": "of analog input 1",
        "zai2": "of analog input 2",
    }

    for z in ["zai1", "zai2"]:
        if z not in ds:
            continue
        ds[z].attrs["standard_name"] = "height"
        ds[z].attrs["units"] = "m"
        ds[z].attrs["positive"] = "up"
        ds[z].attrs["axis"] = "Z"
        ds[z].attrs["long_name"] = f"{long_name} {lnshim[z]}"
        if geopotential_datum_name:
            ds[z].attrs["geopotential_datum_name"] = geopotential_datum_name

    return ds


def ds_drop(ds):
    """
    Drop old DataArrays from Dataset that won't make it into the final .nc file
    """

    todrop = [
        "VEL1",
        "VEL2",
        "VEL3",
        "AMP1",
        "AMP2",
        "AMP3",
        "SNR1",
        "SNR2",
        "SNR3",
        "COR1",
        "COR2",
        "COR3",
        "AnalogInput1",
        "AnalogInput2",
        "Depth",
        "Checksum",
        "ErrorCode",
        "StatusCode",
        "vel1_1277",
        "vel2_1278",
        "vel3_1279",
        "Burst",
        "Battery",
    ]

    if ("AnalogInput1" in ds.attrs) and (ds.attrs["AnalogInput1"].lower() == "true"):
        todrop.remove("AnalogInput1")

    if ("AnalogInput2" in ds.attrs) and (ds.attrs["AnalogInput2"].lower() == "true"):
        todrop.remove("AnalogInput2")

    # drop vars burst and battery if continuous mode to save file space; keep if in burst mode
    if ds.attrs["sample_mode"] == "BURST":
        keep = ["Burst", "Battery"]
        for var in keep:
            todrop.remove(var)

    return ds.drop([t for t in todrop if t in ds.variables])


def scale_analoginput(ds):
    """convert AnalogInput from counts to volts"""

    for v in ["AnalogInput1", "AnalogInput2"]:
        ds[v] = ds[v] * 5 / 65535
        ds = utils.insert_note(
            ds, v, "Converted from counts to volts: volts=counts*5/65535."
        )

    return ds


def associate_z_coord(ds):
    """Associate the appropriate z coordinate to data variables.
    We do this because there are multiple relevant elevations per deployment
    (e.g., velocity and pressure were collected at different elevations,
    and we need to indicate this)"""

    for v in [
        "u_1205",
        "v_1206",
        "w_1204",
        "AGC1_1221",
        "AGC2_1222",
        "AGC3_1223",
        "SNR1",
        "SNR2",
        "SNR3",
        "cor1_1285",
        "cor2_1286",
        "cor3_1287",
    ]:
        if v in ds:
            # pass axis=-1 to add z dim to end for CF compliance
            ds[v] = ds[v].expand_dims("zvel", axis=-1)

    for v in ["P_1ac", "P_1"]:
        if v in ds:
            ds[v] = ds[v].expand_dims("zpres", axis=-1)

    if "AnalogInput1" in ds:
        ds["AnalogInput1"] = ds["AnalogInput1"].expand_dims("zai1", axis=-1)

    if "AnalogInput2" in ds:
        ds["AnalogInput2"] = ds["AnalogInput2"].expand_dims("zai2", axis=-1)

    return ds


def dist_to_boundary(ds):
    """Create range to boundary variable from start/end values"""
    ds["brange"] = (ds["DistProbeStartAvg"] + ds["DistProbeEndAvg"]) / 2
    ds["vrange"] = (ds["DistSVolStartAvg"] + ds["DistSVolEndAvg"]) / 2

    for v in [
        "DistProbeStartAvg",
        "DistProbeEndAvg",
        "DistSVolStartAvg",
        "DistSVolEndAvg",
    ]:
        ds = ds.drop(v)

    return ds


def reshape(ds):

    # find times of first sample in each burst
    t = ds["time"][ds["sample"].values == 1]
    # get corresponding burst number
    b_num = ds["burst"][ds["sample"].values == 1]

    s = []
    t3 = []

    ds.time.encoding.pop("units")

    for i in tqdm(np.arange(0, len(t))):
        t2, samp = np.meshgrid(
            t[i],
            ds["sample"].sel(
                time=slice(
                    t[i],
                    ds.time[ds["burst"] == b_num[i]].max(),
                )
            ),
        )

        s.append(np.array(samp.transpose().flatten()))
        t3.append(np.array(t2.transpose().flatten()))

    s = np.hstack(s)
    t3 = np.hstack(t3)

    ind = pd.MultiIndex.from_arrays((t3, s), names=("new_time", "new_sample"))

    ds = ds.sel(time=slice(t[0], ds["time"][-1])).assign(time=ind).unstack("time")

    ds = ds.drop("sample").rename({"new_time": "time", "new_sample": "sample"})

    return ds


def combine_vars(ds):
    """Combines beam variables (cor, amp, snr) and raw (beam or xyz) velocities into matrix to limit number of variables"""

    ds["amp_avg"] = (ds["AMP1"] + ds["AMP2"] + ds["AMP3"]) / 3

    veldim = ds.attrs["VECCoordinateSystem"].lower()

    if veldim == "xyz":
        veldim = "inst"

    ds["vel"] = xr.DataArray([ds.VEL1, ds.VEL2, ds.VEL3], dims=[veldim, "time"])
    ds["vel"].attrs.update(
        {
            "units": "m s-1",
            "long_name": f"Velocity, {veldim} coordinate system",
        }
    )

    ds["cor"] = xr.DataArray([ds.COR1, ds.COR2, ds.COR3], dims=["beam", "time"]).astype(
        "float32"
    )
    ds["amp"] = xr.DataArray([ds.AMP1, ds.AMP2, ds.AMP3], dims=["beam", "time"]).astype(
        "float32"
    )
    ds["snr"] = xr.DataArray([ds.SNR1, ds.SNR2, ds.SNR3], dims=["beam", "time"])

    return ds


def fill_snr(ds):
    """
    Fill velocity data with corresponding beam snr value threshold
    """
    if "snr_threshold" in ds.attrs:
        Ptxt = str(ds.attrs["snr_threshold"])

        for v in [
            "u_1205",
            "v_1206",
            "w_1204",
        ]:
            if v in ds:
                for bm in ds["beam"].values:
                    ds[v] = ds[v].where(ds.snr.sel(beam=bm) > ds.attrs["snr_threshold"])

                histtext = "Filled velocity data using snr threshold of {} for corresponding beam(s).".format(
                    Ptxt
                )
                ds = utils.insert_note(ds, v, histtext)

        ds = utils.insert_history(ds, histtext)

    return ds


def fill_cor(ds):
    """
    Fill velocity data with corresponding beam correlation value threshold
    """
    if "cor_threshold" in ds.attrs:
        Ptxt = str(ds.attrs["cor_threshold"])

        for v in [
            "u_1205",
            "v_1206",
            "w_1204",
        ]:
            if v in ds:
                for bm in ds["beam"].values:
                    ds[v] = ds[v].where(ds.cor.sel(beam=bm) > ds.attrs["cor_threshold"])

                histtext = "Filled velocity data using cor threshold of {} for corresponding beam(s).".format(
                    Ptxt
                )
                ds = utils.insert_note(ds, v, histtext)

            ds = utils.insert_history(ds, histtext)

    return ds


def remove_vec_globatts(ds):
    """remove unnecessary global attributes from raw instrument file"""

    names = [
        "VECTiltSensor",
        "VECAnalogPowerOutput",
        "VECRecorderSize",
        "VECDeploymentName",
        "VECDNumberOfBeams",
        "VECOutputFormat",
        "VECOutputSync",
        "VECSamplingRate",
        "VECSerialNumber",
    ]

    for att in names:
        del ds.attrs[att]

    return ds


def time_encoding(ds):
    """ensure we don't set dtypes uint for CF compliance"""
    if "units" in ds["time"].encoding:
        ds["time"].encoding.pop("units")

    if ds.attrs["sample_mode"] == "CONTINUOUS":
        if utils.check_time_fits_in_int32(ds, "time"):
            ds["time"].encoding["dtype"] = "i4"
        else:
            print("time variable will not fit in int32; casting to double")
            ds["time"].encoding["dtype"] = "double"
    else:
        if utils.check_time_fits_in_int32(ds, "time"):
            ds["time"].encoding["dtype"] = "i4"

    ds["time"].attrs.update(
        {"standard_name": "time", "axis": "T", "long_name": "time (UTC)"}
    )

    return ds


def var_encoding(ds):
    for var in ds.data_vars:
        if ds[var].dtype == "float64":
            if "dtype" in ds[var].encoding:
                ds[var].encoding.pop("dtype")

            ds[var] = ds[var].astype("float32", copy=False)
            ds[var].encoding["dtype"] = "float32"

    for var in ["burst", "orientation", "sample", "beam"]:
        if var in ds:
            if "dtype" in ds[var].encoding:
                ds[var].encoding.pop("dtype")

            ds[var] = ds[var].astype("int32", copy=False)
            ds[var].encoding["dtype"] = "int32"

    return ds


def create_nc(ds):

    if "prefix" in ds.attrs:
        nc_filename = ds.attrs["prefix"] + ds.attrs["filename"] + "b.nc"
    else:
        nc_filename = ds.attrs["filename"] + "b.nc"

    delayed_obj = ds.to_netcdf(nc_filename, compute=False)

    with ProgressBar():
        delayed_obj.compute()
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])

    print(f"Finished writing data to {nc_filename}")

    return ds


def fill_time_gaps(ds):
    """Fill any gaps in time index to make time even for CONTINUOUS data set before reshaping"""

    if "sample_mode" in ds.attrs:

        if ds.attrs["sample_mode"].upper() == "CONTINUOUS":

            print("Checking for time gaps in CONTINUOUS sampled data")
            sr = ds.attrs["sample_rate"]
            sims = 1 / sr * 1000000000
            pds = (
                int(
                    (ds["time"][-1].values - ds["time"][0].values)
                    / (sims * np.timedelta64(1, "ns"))
                )
                + 1
            )
            idx = pd.date_range(
                str(ds["time"][0].values), periods=pds, freq=f"{sims}ns"
            )

            if len(idx) > len(ds["time"]):
                print(
                    f"Gaps in time index- reindexing to make time index even in CONTINUOUS sampled data using {sims/2} milliseconds tolerance for data variables"
                )

                # make sure time index is unique
                ds = ds.drop_duplicates(dim="time")

                ds = ds.reindex(time=idx, method="nearest", tolerance=f"{sims/2} ms")

        else:
            print(
                f"Not checking for time index gaps for sample mode {ds.attrs['sample_mode']}"
            )

    else:
        print("Not checking for time gaps, sample_mode attribute is missing")

    return ds


def ds_make_average_shape_mi(ds):

    # average_interval is [sec] interval for wave statistics for continuous data
    ds.attrs["average_samples_per_burst"] = int(
        ds.attrs["average_interval"] / ds.attrs["sample_interval"]
    )
    nsamps = ds.attrs["average_samples_per_burst"]

    if len(ds.time) % nsamps != 0:
        ds = ds.isel(time=slice(0, -(int(len(ds.time) % (nsamps)))))

    # create samples & new_time
    x = np.arange(nsamps)
    y = ds["time"][0:-1:nsamps]

    # make new arrays for multi-index
    samp, t = np.meshgrid(x, y)
    s = samp.flatten()
    t3 = t.flatten()

    print(t3.dtype)

    # create multi-index
    ind = pd.MultiIndex.from_arrays((t3, s), names=("new_time", "new_sample"))
    # unstack to make burst shape dataset
    ds = ds.assign_coords(
        xr.Coordinates.from_pandas_multiindex(ind, dim="time")
    ).unstack("time")
    # dsa = dsa.unstack()
    ds = ds.drop_vars("sample").rename({"new_time": "time", "new_sample": "sample"})

    return ds


def drop_more_vars(ds):
    """
    Drop DataArrays from Dataset that are not needed for averaged file
    """

    todrop = ["TransMatrix", "orientmat", "burst", "orientation"]

    return ds.drop([t for t in todrop if t in ds.variables])


def drop_dims(ds):
    """
    Drop dims from Dataset that are not needed for averaged file
    """

    todrop = [
        "inst",
        "earth",
    ]

    return ds.drop([t for t in todrop if t in ds.variables])


def reorder_dims(ds):
    """reorder dimensions for CF compliance"""

    for var in ds.data_vars:
        if "beam" in ds[var].dims:
            ds[var] = ds[var].transpose("beam", "time")

    return ds


def ds_remove_inst_attrs(ds, inst_type="VEC"):
    """Drop instrument attributes from raw Nortek files"""
    rm = []  # initialize list of attrs to be removed

    for j in ds.attrs:
        if re.match(f"^{inst_type}*", j):
            rm.append(j)
    for k in rm:
        del ds.attrs[k]

    return ds


def create_average_nc(ds):

    if "prefix" in ds.attrs:
        nc_filename = ds.attrs["prefix"] + ds.attrs["filename"] + "b-a.nc"
    else:
        nc_filename = ds.attrs["filename"] + "b-a.nc"

    delayed_obj = ds.to_netcdf(nc_filename, compute=False)

    with ProgressBar():
        delayed_obj.compute()

    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])

    print(f"Finished writing data to {nc_filename}")
