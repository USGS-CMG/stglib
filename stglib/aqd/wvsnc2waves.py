import xarray as xr

from ..core import qaqc, utils, waves
from ..vec import nc2waves


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
            if ds.attrs["orientation"] != "UP":
                raise NotImplementedError(
                    "PUV currently only works on UP oriented instruments"
                )

    if dopuv:
        ds = nc2waves.do_puv(ds)

    # ds = utils.create_water_depth(ds)

    # Remove old variables as we just want to keep the wave statistics
    keys = [
        "P_1",
        "P_1ac",
        "sample",
        "Tx_1211",
        "vel1_1277",
        "vel2_1278",
        "vel3_1279",
        "U",
        "V",
        "W",
        "u_1205",
        "v_1206",
        "w_1204",
        "avgamp1",
        "avgamp2",
        "avgamp3",
        "AGC1_1221",
        "AGC2_1222",
        "AGC3_1223",
        "TransMatrix",
        "nrecs",
        "burst",
        "soundspeed",
        "Battery",
        "Hdg_1215",
        "Ptch_1216",
        "Roll_1217",
    ]

    if dopuv:
        keys.remove("sample")

    for k in keys:
        if k in ds:
            ds = ds.drop(k)

    ds = qaqc.drop_vars(ds)

    ds = utils.trim_max_wp(ds)

    ds = utils.trim_min_wh(ds)

    ds = utils.trim_max_wh(ds)

    ds = utils.trim_wp_ratio(ds)

    if dopuv:
        ds = waves.puv_qaqc(ds)

    # Add attrs
    ds = utils.ds_add_wave_attrs(ds)

    nc_filename = ds.attrs["filename"] + "wvs-a.nc"

    ds.to_netcdf(nc_filename, unlimited_dims=["time"])
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])

    print("Done creating", nc_filename)

    return ds
