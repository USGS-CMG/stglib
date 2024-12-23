import re

import numpy as np
import xarray as xr

from ..core import qaqc, utils, waves
from . import nc2waves


def nc_to_diwasp(nc_filename):
    """
    Process burst data to make DIWASP derived wave statistics and spectra
    """
    ds = xr.open_dataset(nc_filename)

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
        ds = nc2waves.make_wave_bursts(ds)

    if "diwasp" not in ds.attrs:
        # if noyt specified choose 'suv' method for signature
        ds.attrs["diwasp"] = "suv"

    if ds.attrs["diwasp"] == "puv" or ds.attrs["diwasp"].lower() == "suv":
        print(f"Running DIWASP using {ds.attrs['diwasp']} input data")
        diwasp = waves.make_diwasp_puv_suv(ds, inst_type="SIG")
        ds = utils.ds_add_pydiwasp_history(ds)
    else:
        raise ValueError(
            f"DIWASP input type {ds.attrs['diwasp']} is not currently supported for {ds.attrs['instument_type']} in stglib"
        )

    for k in [
        "diwasp_frequency",
        "diwasp_direction",
        "diwasp_Hs",
        "diwasp_Tp",
        "diwasp_Tm",
        "diwasp_DTp",
        "diwasp_Dp",
        "diwasp_Dm",
        "diwasp_dspec",
    ]:

        ds[k] = diwasp[k]

    # add diwasp attrs
    for k in diwasp.attrs:
        ds.attrs[k] = diwasp.attrs[k]

    ds = utils.create_water_depth_var(ds)

    for k in [
        "burst",
        "sample",
        "P_1",
        "P_1ac",
        "u_1205",
        "v_1206",
        "w_1204",
        "Tx_1211",
        "SV_80",
        "vel",
        "amp",
        "cor",
        "brangeAST",
        "TransMatrix",
    ]:
        if k in ds:
            ds = ds.drop_vars(k)

    # ds = qaqc.drop_vars(ds)

    ds = nc2waves.drop_unused_dims(ds)

    ds = utils.trim_max_wp(ds)

    ds = utils.trim_min_wh(ds)

    ds = utils.trim_max_wh(ds)

    ds = utils.trim_wp_ratio(ds)

    # Add attrs
    ds = utils.ds_add_wave_attrs(ds)

    # clean-up attrs
    ds = nc2waves.drop_attrs(ds)

    # assign min/max (need to do this after trimming):
    ds = utils.add_min_max(ds)

    # round time to minutes to make sure fits in dtype i4. will be fine for wave burst start times
    ds["time"] = ds["time"].dt.round("min")

    if dodiwasp:
        nc_filename = ds.attrs["filename"] + "s_diwasp-a.nc"

        ds.to_netcdf(nc_filename, unlimited_dims=["time"])
        utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])

        print("Done writing netCDF file", nc_filename)

    return ds
