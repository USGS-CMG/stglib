import numpy as np
import xarray as xr

from ..core import utils, waves
from . import aqdutils


def nc_to_diwasp(nc_filename):
    """
    Process burst data to wave statistics by using DIWASP output from Matlab
    """

    ds = xr.load_dataset(nc_filename)

    mat = xr.load_dataset(ds.attrs["filename"] + "wvs-diwasp.nc")

    ds["frequency"] = xr.DataArray(mat["frequency"], dims=("frequency"))

    # Note we convert the polar/cartesian coordinates provided by DIWASP to
    # compass coordinates
    dirs = waves.polar2compass(waves.to2from(mat["direction"]))

    # Return only unique directions (diwasp has double directions...)
    # FIXME: Why is this?
    _, idx = np.unique(dirs, return_index=True)

    ds["direction"] = xr.DataArray(dirs[idx], dims=("direction"))

    ds["dspec"] = xr.DataArray(
        mat["dspec"][:, idx, :], dims=("time", "direction", "frequency")
    )

    pspec = np.trapz(ds["dspec"].values, x=ds["direction"], axis=1)

    m0 = waves.make_moment(ds["frequency"].values, pspec, 0)
    m2 = waves.make_moment(ds["frequency"].values, pspec, 2)

    ds["wh_4061"] = xr.DataArray(waves.make_Hs(m0), dims="time")

    ds["wp_4060"] = xr.DataArray(waves.make_Tm(m0, m2), dims="time")

    for k in ["wp_peak"]:
        ds[k] = xr.DataArray(mat[k], dims="time")

    for k in ["wvdir", "dwvdir"]:
        ds[k] = xr.DataArray(waves.polar2compass(waves.to2from(mat[k])), dims="time")

    ds["pspec"] = xr.DataArray(pspec, dims=("time", "frequency"))

    ds["wd_4062"] = xr.DataArray(
        waves.polar2compass(
            waves.to2from(
                waves.make_mwd(
                    ds["frequency"].values, ds["direction"].values, ds["dspec"].values
                )
            )
        ),
        dims="time",
    )

    # Remove old variables as we just want to keep the wave statistics
    for v in [
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
        "avgamp1",
        "avgamp2",
        "avgamp3",
        "AGC1_1221",
        "AGC2_1222",
        "AGC3_1223",
        "TransMatrix",
        "soundspeed",
        "Battery",
        "Hdg_1215",
        "Ptch_1216",
        "Roll_1217",
        "Bat_106",
        "Depth",
    ]:
        if v in ds:
            ds = ds.drop(v)

    # remove depth from bin_depth
    ds["bin_depth"] = ds["bin_depth"].squeeze(dim="depth")

    ds = utils.trim_max_wp(ds)

    ds = utils.trim_min_wh(ds)

    ds = utils.trim_max_wh(ds)

    ds = utils.trim_wp_ratio(ds)

    # Add attrs
    ds = aqdutils.ds_add_attrs(ds, waves=True)

    ds = utils.ds_add_diwasp_history(ds)

    nc_filename = ds.attrs["filename"] + "wvs_diwasp-cal.nc"

    print("Writing to netCDF")

    ds.to_netcdf(nc_filename)
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])

    print("Done creating", nc_filename)

    return ds
