import numpy as np
import xarray as xr

from ..core import qaqc, utils, waves


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

    if dopuv:
        ds = do_puv(ds)

    # keep burst mean P_1 and P_1ac for reference
    for k in ["P_1", "P_1ac"]:
        if k in ds:
            ds[k] = ds[k].mean(dim="sample", keep_attrs=True)

    for k in [
        "burst",
        "brange",
        "sample",
        "u_1205",
        "v_1206",
        "w_1204",
        "vel1_1277",
        "vel2_1278",
        "vel3_1279",
        "vel",
        "vrange",
        "amp",
        "AGC_1202",
        "AGC1_1221",
        "AGC2_1222",
        "AGC3_1223",
        "SNR1",
        "SNR2",
        "SNR3",
        "snr",
        "cor1_1285",
        "cor2_1286",
        "cor3_1287",
        "cor",
        "AnalogInput1",
        "AnalogInput2",
        "Hdg_1215",
        "Ptch_1216",
        "Roll_1217",
        "Bat_106",
        "Tx_1211",
        "SV_80",
        "orientation",
        "orientmat",
    ]:
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

    # assign min/max (need to do this after trimming):
    ds = utils.add_min_max(ds)

    nc_filename = ds.attrs["filename"] + "s-a.nc"

    ds.to_netcdf(nc_filename, unlimited_dims=["time"])
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])
    print("Done writing netCDF file", nc_filename)

    return ds


def do_puv(ds):
    print("Running puv_quick")

    for k in ["pressure_sensor_height", "velocity_sample_volume_height"]:
        if k not in ds.attrs:
            raise KeyError(f"{k} must be specified to run PUV")

    N, M = np.shape(ds["u_1205"].squeeze())

    if "puv_first_frequency_cutoff" in ds.attrs:
        first_frequency_cutoff = ds.attrs["puv_first_frequency_cutoff"]
    else:
        first_frequency_cutoff = 1 / 10

    if "puv_last_frequency_cutoff" in ds.attrs:
        last_frequency_cutoff = ds.attrs["puv_last_frequency_cutoff"]
    else:
        last_frequency_cutoff = 1 / 2.5

    if "P_1ac" in ds:
        pvar = "P_1ac"
    else:
        pvar = "P_1"

    u = ds["u_1205"].squeeze()
    v = ds["v_1206"].squeeze()

    psh = ds.attrs["pressure_sensor_height"]
    huv = ds.attrs["velocity_sample_volume_height"]

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
