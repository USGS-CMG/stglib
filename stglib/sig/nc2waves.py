import re
import warnings

import numpy as np
import pandas as pd
import xarray as xr

from ..core import utils, waves
from .cdf2nc import add_water_level


def nc_to_waves(nc_filename, salwtemp=None):
    """
    Process burst data to wave statistics
    """
    ds = xr.open_dataset(nc_filename, chunks={"time": 200000, "z": 5})

    # check to see if need to make wave burst from continuous data
    if (ds.attrs["sample_mode"] == "CONTINUOUS") and ("wave_interval" in ds.attrs):
        # check for wave_start_time attrs
        if "wave_start_time" in ds.attrs:
            print(
                f"trimming continuous data to start at users specified wave start time {ds.attrs['wave_start_time']}"
            )
            ds = ds.sel(
                time=slice(np.datetime64(ds.attrs["wave_start_time"]), ds["time"][-1])
            )
        # make wave burst ncfile from continuous data if wave_interval is specified
        ds = make_wave_bursts_mi(ds)

    if "BURST" in ds.attrs["sample_mode"]:
        if "wave_interval" not in ds.attrs:
            if "burst_interval" in ds.attrs:
                ds.attrs["wave_interval"] = ds.attrs["burst_interval"]
            else:
                ds.attrs["wave_interval"] = int(
                    (ds["time"][1] - ds["time"][0]).dt.round("s").values
                    / np.timedelta64(1, "s")
                )

    # check to see if wave duration is specified and if so trim burst samples accordingly
    if "wave_duration" in ds.attrs:
        if ds.attrs["wave_duration"] <= ds.attrs["wave_interval"]:
            nsamps = int(ds.attrs["wave_duration"] * ds.attrs["sample_rate"])
            ds = ds.isel(sample=slice(0, nsamps))
            ds.attrs["wave_samples_per_burst"] = nsamps

        else:
            warnings.warn(
                f"wave_duration {ds.attrs.wave_duration} cannot be greater than wave_interval {ds.attrs.wave_interval}, this configuration setting will be ignored"
            )

    spec = waves.make_waves_ds(ds)
    ds = utils.ds_add_waves_history(ds)

    for k in ["wp_peak", "wh_4061", "wp_4060", "pspec"]:
        ds[k] = spec[k]

    dopuv = False
    if "puv" in ds.attrs:
        if ds.attrs["puv"].lower() == "true":
            dopuv = True

    if dopuv:
        ds = do_puv(ds)

    if "pressure_sensor_height" in ds.attrs:
        ds = utils.create_water_depth_var(
            ds.mean(dim="sample", keep_attrs=True),
            psh="pressure_sensor_height",
            salwtemp=salwtemp,
        )
    else:
        ds = utils.create_water_depth_var(
            ds.mean(dim="sample", keep_attrs=True), salwtemp=salwtemp
        )

    ds = add_water_level(ds, salwtemp=salwtemp)

    for k in [
        "burst",
        "sample",
        "P_1",
        "P_1ac",
        "u_1205",
        "v_1206",
        "w_1204",
        "w2_1204",
        "Tx_1211",
        "SV_80",
        "vel",
        "amp",
        "cor",
        "brangeAST",
        "TransMatrix",
        "ast_quality",
    ]:

        if k in ds:
            ds = ds.drop_vars(k)

    # ds = qaqc.drop_vars(ds)

    ds = make_waves_vdims(ds, wtype="nc2waves")

    ds = drop_unused_dims(ds)

    ds = utils.trim_max_wp(ds)

    ds = utils.trim_min_wh(ds)

    ds = utils.trim_max_wh(ds)

    ds = utils.trim_wp_ratio(ds)

    if dopuv:
        ds = waves.puv_qaqc(ds)

    # Add attrs
    ds = utils.ds_add_wave_attrs(ds)

    # clean-up attrs
    ds = drop_attrs(ds)

    # assign min/max (need to do this after trimming):
    ds = utils.add_min_max(ds)

    ds = utils.add_delta_t(ds)

    # round time to minutes to make sure fits in dtype i4. will be fine for wave burst start times
    ds["time"] = ds["time"].dt.round("min")

    if dopuv:
        nc_filename = ds.attrs["filename"] + "s_puvq-a.nc"
    else:
        nc_filename = ds.attrs["filename"] + "s-a.nc"

    ds.to_netcdf(nc_filename, unlimited_dims=["time"])
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])

    print("Done writing netCDF file", nc_filename)

    return ds


def nc_to_diwasp(nc_filename, salwtemp=None):
    """
    Process burst data to make DIWASP derived wave statistics and spectra
    """

    ds = xr.open_dataset(nc_filename, chunks={"time": 200000, "z": 5})

    # check to see if need to make wave burst from continuous data
    if (ds.attrs["sample_mode"] == "CONTINUOUS") and ("wave_interval" in ds.attrs):
        # check for wave_start_time attrs
        if "wave_start_time" in ds.attrs:
            print(
                f"trimming continuous data to start at users specified wave start time {ds.attrs['wave_start_time']}"
            )
            ds = ds.sel(
                time=slice(np.datetime64(ds.attrs["wave_start_time"]), ds["time"][-1])
            )
        # make wave burst ncfile from continuous data if wave_interval is specified
        ds = make_wave_bursts_mi(ds)

    if "BURST" in ds.attrs["sample_mode"]:
        if "wave_interval" not in ds.attrs:
            if "burst_interval" in ds.attrs:
                ds.attrs["wave_interval"] = ds.attrs["burst_interval"]
            else:
                ds.attrs["wave_interval"] = int(
                    (ds["time"][1] - ds["time"][0]).dt.round("s").values
                    / np.timedelta64(1, "s")
                )

    # check to see if wave duration is specified and if so trim burst samples accordingly
    if "wave_duration" in ds.attrs:
        if ds.attrs["wave_duration"] <= ds.attrs["wave_interval"]:
            nsamps = int(ds.attrs["wave_duration"] * ds.attrs["sample_rate"])
            ds = ds.isel(sample=slice(0, nsamps))
            ds.attrs["wave_samples_per_burst"] = nsamps

        else:
            warnings.warn(
                f"wave_duration {ds.attrs.wave_duration} cannot be greater than wave_interval {ds.attrs.wave_interval}, this configuration setting will be ignored"
            )

    if "diwasp" not in ds.attrs:
        # if not specified choose 'suv' method for signature
        ds.attrs["diwasp"] = "suv"
        print("diwasp processing type was not specified, defaulting to 'suv' method")

    if "diwasp" in ds.attrs:
        data_type = ds.attrs["diwasp"]
        if data_type not in ["suv", "puv", "elev", "pres", "optimized", "optimized-nd"]:
            raise ValueError(
                f"data type {data_type} is not recognized current options are ['suv', 'puv', 'elev', 'pres', 'optimized', 'optimized-nd']"
            )

    if "diwasp_ibin" in ds.attrs:
        ibin = ds.attrs["diwasp_ibin"]
    else:
        ibin = 0

    if data_type == "puv" or data_type == "suv" or data_type == "optimized":
        print(f"Running DIWASP using {data_type} input data")
        layout = make_diwasp_layout(ds, data_type=data_type, ibin=ibin)
        diwasp = waves.make_diwasp_puv_suv(
            ds, data_type=data_type, layout=layout, ibin=ibin
        )
        ds = utils.ds_add_pydiwasp_history(ds)

    elif data_type == "elev" or data_type == "pres" or data_type == "optimized-nd":
        print(f"Running DIWASP using {ds.attrs['diwasp']} input data")
        layout = make_diwasp_layout(ds, data_type=data_type)
        diwasp = waves.make_diwasp_elev_pres(ds, data_type=data_type, layout=layout)
        ds = utils.ds_add_pydiwasp_history(ds)
    else:
        raise ValueError(
            f"DIWASP input type {ds.attrs['diwasp']} is not currently supported for {ds.attrs['instument_type']} in stglib"
        )

    for k in [
        "diwasp_frequency",
        "diwasp_direction",
        "diwasp_hs",
        "diwasp_tp",
        "diwasp_tm",
        "diwasp_dtp",
        "diwasp_dp",
        "diwasp_dm",
        "diwasp_dspec",
        "diwasp_fspec",
        "diwasp_type",
    ]:

        if k in diwasp:
            ds[k] = diwasp[k]

    # add diwasp attrs
    for k in diwasp.attrs:
        ds.attrs[k] = diwasp.attrs[k]

    # rename Fspec based on input datatype
    if "optimized" not in data_type:
        ds = utils.rename_diwasp_fspec(ds)

    if "pressure_sensor_height" in ds.attrs:
        ds = utils.create_water_depth_var(
            ds.mean(dim="sample", keep_attrs=True),
            psh="pressure_sensor_height",
            salwtemp=salwtemp,
        )
    else:
        ds = utils.create_water_depth_var(
            ds.mean(dim="sample", keep_attrs=True), salwtemp=salwtemp
        )

    ds = add_water_level(ds, salwtemp=salwtemp)

    for k in [
        "burst",
        "sample",
        "P_1",
        "P_1ac",
        "u_1205",
        "v_1206",
        "w_1204",
        "w2_1204",
        "Tx_1211",
        "SV_80",
        "vel",
        "amp",
        "cor",
        "brangeAST",
        "TransMatrix",
        "ast_quality",
    ]:

        if k in ds:
            ds = ds.drop_vars(k)

    # ds = qaqc.drop_vars(ds)

    # use "epic" names
    ds.attrs["diwasp_names"] = "epic"
    ds = utils.rename_diwasp_wave_vars(ds)

    ds = make_waves_vdims(ds, wtype="nc2diwasp")

    ds = drop_unused_dims(ds)

    ds = utils.trim_max_wp(ds)

    ds = utils.trim_min_wh(ds)

    ds = utils.trim_max_wh(ds)

    ds = utils.trim_wp_ratio(ds)

    # Add attrs
    ds = utils.ds_add_wave_attrs(ds)

    # clean-up attrs
    ds = drop_attrs(ds)

    # assign min/max (need to do this after trimming):
    ds = utils.add_min_max(ds)

    ds = utils.add_delta_t(ds)

    # round time to minutes to make sure fits in dtype i4. will be fine for wave burst start times
    ds["time"] = ds["time"].dt.round("min")

    nc_filename = ds.attrs["filename"] + "s_diwasp-a.nc"

    ds.to_netcdf(nc_filename, unlimited_dims=["time"])
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])

    print("Done writing netCDF file", nc_filename)

    return ds


def make_diwasp_layout(ds, data_type=None, ibin=None):
    """
    Make layout for DIWASP wave processing input
    """

    if data_type == "suv" or data_type == "puv" or data_type == "optimized":
        # datatypes = ["pres", "velx", "vely"]
        sxyz = [0, 0, ds.attrs["initial_instrument_height"]]
        if ds.attrs["orientation"].lower() == "up":
            velz = ds.attrs["initial_instrument_height"] + ds["bindist"][ibin].values

        elif ds.attrs["orientation"].lower() == "down":
            velz = ds.attrs["initial_instrument_height"] - ds["bindist"][ibin].values

        uxyz = vxyz = [0, 0, velz]

        layout = np.array([sxyz, uxyz, vxyz]).T

    elif data_type == "elev" or data_type == "pres":
        # datatypes=['elev']
        sxyz = np.atleast_2d([0, 0, ds.attrs["initial_instrument_height"]]).T
        layout = sxyz

    return layout


def make_wave_bursts_mi(
    ds,
    wave_vars=[
        "sample",
        "P_1",
        "P_1ac",
        "Tx_1211",
        "brangeAST",
        "ast_quality",
        "u_1205",
        "v_1206",
        "w_1204",
        "vel",
        "cor",
        "amp",
    ],
):
    """
    Reshape CONTINUOUS data set into burst shape using multi-indexing for purpose of wave analysis
    """
    for k in ds.data_vars:
        if k not in wave_vars:
            ds = ds.drop_vars(k)

    # average_interval is [sec] interval for wave statistics for continuous data
    ds.attrs["wave_samples_per_burst"] = int(
        ds.attrs["wave_interval"] / ds.attrs["sample_interval"]
    )
    nsamps = ds.attrs["wave_samples_per_burst"]

    if len(ds.time) % nsamps != 0:
        ds = ds.isel(time=slice(0, -(int(len(ds.time) % (nsamps)))))

    # create samples & new_time
    x = np.arange(nsamps)
    y = ds["time"][0:-1:nsamps]

    # make new arrays for multi-index
    samp, t = np.meshgrid(x, y)
    s = samp.flatten()
    t3 = t.flatten()

    # create multi-index
    ind = pd.MultiIndex.from_arrays((t3, s), names=("new_time", "new_sample"))
    # unstack to make burst shape dataset
    ds = ds.assign_coords(
        xr.Coordinates.from_pandas_multiindex(ind, dim="time")
    ).unstack("time")
    # dsa = dsa.unstack()
    if "sample" in ds.data_vars:
        ds = ds.drop_vars("sample").rename({"new_time": "time", "new_sample": "sample"})
    else:
        ds = ds.rename({"new_time": "time", "new_sample": "sample"})

    return ds


def make_wave_bursts(ds):
    # wave_interval is [sec] interval for wave statistics for continuous data
    ds.attrs["wave_samples_per_burst"] = int(
        ds.attrs["wave_interval"] / ds.attrs["sample_interval"]
    )

    r = np.shape(ds.P_1)[0]
    mod = r % ds.attrs["wave_samples_per_burst"]
    if mod:
        print(
            "Number of rows is not a multiple of samples_per_burst; truncating to last full burst"
        )
        ds = ds.sel(time=ds.time[0:-mod])

    ds.time.encoding.pop("dtype")

    ds["timenew"] = xr.DataArray(
        ds.time[0 :: int(ds.attrs["wave_samples_per_burst"])].values, dims="timenew"
    )

    ds["samplenew"] = xr.DataArray(
        range(ds.attrs["wave_samples_per_burst"]), dims="samplenew"
    )

    for v in ["P_1", "P_1ac", "Tx_1211", "brangeAST", "ast_quality"]:
        if v in ds:
            attrsbak = ds[v].attrs
            ds[v] = xr.DataArray(
                np.reshape(ds[v].values, (-1, int(ds.attrs["wave_samples_per_burst"]))),
                dims=["timenew", "samplenew"],
            )
            ds[v].attrs = attrsbak

    for v in ["u_1205", "v_1206", "w_1204"]:
        if v in ds:
            attrsbak = ds[v].attrs
            # check dims
            if ds[v].dims == ("time", "z"):
                ds[v] = xr.DataArray(
                    np.reshape(
                        ds[v].values,
                        (-1, int(ds.attrs["wave_samples_per_burst"]), len(ds["z"])),
                    ).transpose(0, 2, 1),
                    dims=["timenew", "z", "samplenew"],
                )
                # u=np.reshape(ds['u_1205'].values, (-1,int(ds.attrs["samples_per_burst"]), len(ds['z']))).transpose(0,2,1)
                ds[v].attrs = attrsbak
            else:
                raise ValueError(
                    f"{v} dimensions {ds[v].dims} are not as required ('time','z') to shape into wave burst"
                )

    """
    for v in ["vel", "cor", "amp"]:
        if v in ds:
            attrsbak = ds[v].attrs
            if ds[v].dims == ("beam", "time", "z"):
                ds[v] = xr.DataArray(
                    np.reshape(
                        ds[v].values,
                        (
                            len(ds["beam"]),
                            -1,
                            int(ds.attrs["wave_samples_per_burst"]),
                            len(ds["z"]),
                        ),
                    ).transpose(0, 1, 3, 2),
                    dims=["beam", "timenew", "z", "samplenew"],
                )
                # u=np.reshape(ds['u_1205'].values, (-1,int(ds.attrs["samples_per_burst"]), len(ds['z']))).transpose(0,2,1)
                ds[v].attrs = attrsbak
            else:
                raise ValueError(
                    f"{v} dimensions {ds[v].dims} are not as required ('beam','time','z') to shape into wave burst"
                )
    """
    ds = ds.rename({"time": "timeold"})
    ds = ds.rename({"timenew": "time"})
    ds = ds.drop_dims("timeold")

    if "sample" in ds:
        ds = ds.rename({"sample": "sampleold"})
        ds = ds.rename({"samplenew": "sample"})
        ds = ds.drop_vars("sampleold")
    else:
        ds = ds.rename({"samplenew": "sample"})

    return ds


def drop_unused_dims(ds):
    """only keep dims that will be in the final files"""
    thedims = [
        "frequency",
        "direction",
        "latitude",
        "longitude",
        "z",
        "depth",
        "zsen",
        "depthsen",
    ]

    for v in ds.data_vars:
        for x in ds[v].dims:
            thedims.append(x)

    for x in ds.dims:
        if x not in thedims:
            ds = ds.drop_vars(x)

    return ds


def drop_attrs(ds):
    """Drop some global attrs"""

    rm = []  # initialize list of attrs to be removed
    exclude = []  # initialize attrs to exclude

    for j in ds.attrs:
        if re.match("^SIG", j):
            rm.append(j)

    for k in rm:
        if k not in exclude:
            del ds.attrs[k]

    return ds


def do_puv(ds):
    print("Running puv_quick")

    for k in ["initial_instrument_height"]:
        if k not in ds.attrs:
            raise KeyError(f"{k} must be specified to run PUV")

    if "puv_bin" in ds.attrs:
        ibin = ds.attrs["puv_bin"]
    else:
        ibin = 0

    N, M = np.shape(ds["u_1205"].isel(z=ibin).squeeze())

    if "puv_first_frequency_cutoff" in ds.attrs:
        first_frequency_cutoff = ds.attrs["puv_first_frequency_cutoff"]
    else:
        first_frequency_cutoff = 1 / 33

    if "puv_last_frequency_cutoff" in ds.attrs:
        last_frequency_cutoff = ds.attrs["puv_last_frequency_cutoff"]
    else:
        last_frequency_cutoff = 1 / 3.3

    if "P_1ac" in ds:
        pvar = "P_1ac"
    else:
        pvar = "P_1"

    u = ds["u_1205"].isel(z=ibin).squeeze()
    v = ds["v_1206"].isel(z=ibin).squeeze()

    psh = ds.attrs["initial_instrument_height"]
    if ds.attrs["orientation"].lower() == "up":
        huv = ds.attrs["initial_instrument_height"] + ds["bindist"][ibin].values
    elif ds.attrs["orientation"].lower() == "down":
        huv = ds.attrs["initial_instrument_height"] - ds["bindist"][ibin].values

    ds = waves.call_puv_quick_vectorized(
        ds,
        pvar=pvar,
        u=u,
        v=v,
        psh=psh,
        huv=huv,
        first_frequency_cutoff=first_frequency_cutoff,
        last_frequency_cutoff=last_frequency_cutoff,
    )

    return ds


def make_waves_vdims(ds, wtype="nc2waves"):
    """Modify vertical dimensions for waves output"""
    if wtype == "nc2waves":
        if "puv" in ds.attrs and ds.attrs["puv"].lower() == "true":
            data_type = "puv"
            if "puv_bin" in ds.attrs:
                ibin = ds.attrs["puv_bin"]
            else:
                ibin = 0
        else:
            data_type = "pres"

    elif wtype == "nc2diwasp":
        data_type = ds.attrs["diwasp"]

        if data_type in ["suv", "puv", "optimized"]:
            if "diwasp_ibin" in ds.attrs:
                ibin = ds.attrs["diwasp_ibin"]
            else:
                ibin = 0

    attrsz = ds["z"].attrs
    attrsdep = ds["depth"].attrs
    notez = None
    notedep = None
    if data_type in ["suv", "puv", "optimized"]:
        # set vertical dimension for veloctiy measurement used in suv, puv calcs
        if "z" in ds.dims:
            ds["z"] = xr.DataArray([ds["z"][ibin].values], dims="z")
            notez = "height to wave velocity measurement"

        if "depth" in ds.dims:
            ds["depth"] = xr.DataArray([ds["depth"][ibin].values], dims="depth")
            notedep = "depth to wave velocity measurement"

    elif data_type in ["pres", "elev", "optimized-nd"]:
        # swap vertical dimension to sensor location if not doing puv
        bin0dist = ds["bindist"][0].values
        if "zsen" in ds.dims:
            ds["z"] = xr.DataArray(ds["zsen"].values, dims="z")
        else:
            if ds.attrs["orientation"].upper() == "DOWN":
                ds["z"] = xr.DataArray([ds["z"][0].values + bin0dist], dims="z")
            if ds.attrs["orientation"].upper() == "UP":
                ds["z"] = xr.DataArray([ds["z"][0].values - bin0dist], dims="z")

        if "depthsen" in ds.dims:
            ds["depth"] = xr.DataArray(ds["depthsen"].values, dims="depth")

        else:
            if ds.attrs["orientation"].upper() == "DOWN":
                ds["depth"] = xr.DataArray(
                    [ds["depth"][0].values - bin0dist], dims="depth"
                )
            if ds.attrs["orientation"].upper() == "UP":
                ds["depth"] = xr.DataArray(
                    [ds["depth"][0].values + bin0dist], dims="depth"
                )

        vdimsdrop = ["zsen", "depthsen", "bindist"]
        for k in vdimsdrop:
            if k in ds.dims:
                ds = ds.drop_vars(k)

    ds["z"].attrs = attrsz
    ds["depth"].attrs = attrsdep
    if notez:
        ds["z"].attrs["note"] = notez
    if notedep:
        ds["depth"].attrs["note"] = notedep

    return ds
