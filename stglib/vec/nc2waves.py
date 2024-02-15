import numpy as np
import xarray as xr

from ..core import utils, waves


def nc_to_waves(nc_filename):
    """
    Process burst data to wave statistics
    """

    ds = xr.load_dataset(nc_filename)

    spec = waves.make_waves_ds(ds)

    for k in ["wp_peak", "wh_4061", "wp_4060", "pspec"]:
        ds[k] = spec[k]

    dopuv = False
    if "puv" in ds.attrs:
        if ds.attrs["puv"].lower() == "true":
            dopuv = True

    # keep burst mean P_1 and P_1ac for reference
    for k in ["P_1", "P_1ac"]:
        if k in ds:
            ds[k] = ds[k].mean(dim="sample", keep_attrs=True)

    for k in [
        "burst",
        "u_1205",
        "v_1206",
        "w_1204",
        "vel1_1277",
        "vel2_1278",
        "vel3_1279",
        "AGC1_1221",
        "AGC2_1222",
        "AGC3_1223",
        "SNR1",
        "SNR2",
        "SNR3",
        "cor1_1285",
        "cor2_1286",
        "cor3_1287",
        "AnalogInput1",
        "AnalogInput2",
        "Hdg_1215",
        "Ptch_1216",
        "Roll_1217",
        "Bat_106",
        "Tx_1211",
        "SV_80",
    ]:
        if k in ds:
            ds = ds.drop(k)

    ds = utils.trim_max_wp(ds)

    ds = utils.trim_min_wh(ds)

    ds = utils.trim_max_wh(ds)

    ds = utils.trim_wp_ratio(ds)

    if dopuv:
        ds = waves.puv_qaqc(ds)

    # Add attrs
    ds = utils.ds_add_wave_attrs(ds)

    # assign min/max (need to do this after trimming):
    ds = utils.add_min_max(ds)

    nc_filename = ds.attrs["filename"] + "s-a.nc"

    ds.to_netcdf(nc_filename, unlimited_dims=["time"])
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])
    print("Done writing netCDF file", nc_filename)

    return ds
