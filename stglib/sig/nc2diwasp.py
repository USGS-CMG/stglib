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

    if "diwasp" in ds.attrs:
        if "suv" in ds.attrs["diwasp"]:
            data_type = "suv"
        elif "puv" in ds.attrs["diwasp"]:
            data_type = "puv"
        elif "elev" in ds.attrs["diwasp"]:
            data_type = "elev"
        elif "pres" in ds.attrs["diwasp"]:
            data_type = "pres"

    if "diwasp_bin" in ds.attrs:
        ibin = ds.attrs["diwasp_bin"]
    else:
        ibin = 0

    layout = make_layout(ds, data_type=data_type, ibin=ibin)

    if data_type == "puv" or data_type == "suv":
        print(f"Running DIWASP using {data_type} input data")
        layout = make_layout(ds, data_type=data_type, ibin=ibin)
        diwasp = waves.make_diwasp_puv_suv(
            ds, data_type=data_type, layout=layout, ibin=ibin
        )
        ds = utils.ds_add_pydiwasp_history(ds)

    elif ds.attrs["diwasp"] == "elev" or ds.attrs["diwasp"].lower() == "pres":
        print(f"Running DIWASP using {ds.attrs['diwasp']} input data")
        layout = make_layout(ds, data_type=data_type)
        diwasp = waves.make_diwasp_elev_pres(ds, data_type=data_type, layout=layout)
        ds = utils.ds_add_pydiwasp_history(ds)
    else:
        raise ValueError(
            f"DIWASP input type {ds.attrs['diwasp']} is not currently supported for {ds.attrs['instument_type']} in stglib"
        )

    for k in [
        "diwaspFrequency",
        "diwaspDirection",
        "diwaspHs",
        "diwaspTp",
        "diwaspTm",
        "diwaspDTp",
        "diwaspDp",
        "diwaspDm",
        "diwaspDspec",
        "diwaspFspec",
    ]:

        if k in diwasp:
            ds[k] = diwasp[k]

    # add diwasp attrs
    for k in diwasp.attrs:
        ds.attrs[k] = diwasp.attrs[k]

    # rename Fspec based on input datatype
    ds = rename_diwasp_Fspec(ds)

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

    if "diwasp_names" in ds.attrs:
        ds = ds_rename_wave_vars(ds)

    # Add attrs
    ds = utils.ds_add_wave_attrs(ds)

    # clean-up attrs
    ds = nc2waves.drop_attrs(ds)

    # assign min/max (need to do this after trimming):
    ds = utils.add_min_max(ds)

    # round time to minutes to make sure fits in dtype i4. will be fine for wave burst start times
    ds["time"] = ds["time"].dt.round("min")

    nc_filename = ds.attrs["filename"] + "s_diwasp-a.nc"

    ds.to_netcdf(nc_filename, unlimited_dims=["time"])
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])

    print("Done writing netCDF file", nc_filename)

    return ds


def make_layout(ds, data_type=None, ibin=None):
    """
    Make layout for DIWASP wave processing input
    """
    if data_type:

        # get layout
        if data_type == "puv":
            # datatypes = ["pres", "velx", "vely"]
            pxyz = [0, 0, ds.attrs["initial_instrument_height"]]
            if ds.attrs["orientation"].lower() == "up":
                uxyz = [
                    0,
                    0,
                    ds.attrs["initial_instrument_height"] + ds["bindist"][ibin].values,
                ]
            elif ds.attrs["orientation"].lower() == "down":
                uxyz = [
                    0,
                    0,
                    ds.attrs["initial_instrument_height"] - ds["bindist"][ibin].values,
                ]

            vxyz = uxyz

            layout = np.array([pxyz, uxyz, vxyz])

        elif data_type == "suv":
            # datatypes = ["elev", "velx", "vely"]
            sxyz = [0, 0, ds.attrs["initial_instrument_height"]]
            if ds.attrs["orientation"].lower() == "up":
                uxyz = [
                    0,
                    0,
                    ds.attrs["initial_instrument_height"] + ds["bindist"][ibin].values,
                ]
            elif ds.attrs["orientation"].lower() == "down":
                uxyz = [
                    0,
                    0,
                    ds.attrs["initial_instrument_height"] - ds["bindist"][ibin].values,
                ]

            vxyz = uxyz

            layout = np.array([sxyz, uxyz, vxyz])

        elif data_type == "elev":
            # datatypes=['elev']
            sxyz = np.atleast_2d([0, 0, ds.attrs["initial_instrument_height"]]).T
            layout = sxyz

        elif data_type == "pres":
            # datatypes =['pres']
            pxyz = np.atleast_2d([0, 0, ds.attrs["initial_instrument_height"]]).T
            layout = pxyz

        return layout


def ds_rename_wave_vars(ds):
    # rename wave vars to user specified convention
    if ds.attrs["diwasp_names"].lower() == "epic":
        varnames = {
            "diwaspFrequency": "frequency",
            "diwaspDirection": "direction",
            "diwaspHs": "wh_4061",
            "diwaspTp": "wp_peak",
            "diwaspTm": "wp_4060",
            "diwaspDTp": "wvdir",
            "diwaspDp": "dwvdir",
            "diwaspDm": "wd_4062",
            "diwaspPspec": "pspec",
            "diwaspASTspec": "sspec",
            "diwaspVspec": "vspec",
            "diwaspDspec": "dspec",
        }

        # check to make sure they exist before trying to rename
        newvars = {}
        for k in varnames:
            if k in ds:
                newvars[k] = varnames[k]
    else:
        return ds

    return ds.rename(newvars)


def rename_diwasp_Fspec(diwasp):
    # if it exists rename diwaspFspec to first input datatype
    if (
        "diwaspFspec" not in diwasp.data_vars
    ):  # check to make sure Fspec is in data_vars
        return
    else:
        if "diwasp_inputs" in diwasp.attrs:

            inputs = diwasp.attrs["diwasp_inputs"]
            newname = {}
            if inputs[0] == "elev":
                newname = {"diwaspFspec": "diwaspASTspec"}
            elif inputs[0] == "pres":
                newname = {"diwaspFspec": "diwaspPspec"}
            elif (
                inputs[0] == "velx"
                or diwasp.attrs["diwasp_inputs"][0] == "vely"
                or diwasp.attrs["diwasp_inputs"][0] == "velz"
                or diwasp.attrs["diwasp_inputs"][0] == "radial"
            ):
                newname = {"diwaspFspec": "diwaspVspec"}
            else:
                return diwasp
        else:
            return diwasp

    return diwasp.rename(newname)
