from pathlib import Path

import xarray as xr

from ..aqd import aqdutils
from ..core import qaqc, utils


def cdf_to_nc(cdf_filename, atmpres=None):
    """
    Load a raw .cdf file and generate a processed .nc file
    """

    # Load raw .cdf data
    ds = aqdutils.load_cdf(cdf_filename, atmpres=atmpres)

    # Clip data to in/out water times or via good_ens
    ds = utils.clip_ds(ds)

    ds = utils.create_nominal_instrument_depth(ds)

    ds = aqdutils.set_orientation(ds, inst_type="RDI")

    for var, beam in zip(
        ["cor1_1285", "cor2_1286", "cor3_1287", "cor4_1288"], [1, 2, 3, 4]
    ):
        ds[var] = xr.DataArray(ds["corr"].sel(beam=beam), dims=("time", "bindist"))
        ds[var].attrs["long_name"] = (f"Correlation Beam {beam}",)
        ds[var].attrs["units"] = "1"

    for var, beam in zip(
        ["AGC1_1221", "AGC2_1222", "AGC3_1223", "AGC4_1224"], [1, 2, 3, 4]
    ):
        ds[var] = xr.DataArray(ds["int"].sel(beam=beam), dims=("time", "bindist"))
        ds[var].attrs["long_name"] = f"Echo Intensity (AGC) Beam {beam}"
        ds[var].attrs[
            "standard_name"
        ] = "signal_intensity_from_multibeam_acoustic_doppler_velocity_sensor_in_sea_water"
        ds[var].attrs["units"] = "1"

    # Not dealing with percent good for now.
    # for var, beam in zip(
    #     ["PGd1_1241", "PGd2_1242", "PGd3_1243", "PGd4_1244"], [1, 2, 3, 4]
    # ):
    #     ds[var] = xr.DataArray(ds["pg"].sel(beam=beam), dims=("time", "bindist"))
    #     ds[var].attrs["long_name"] = f"Percent Good Beam {beam}"

    for var, beam in zip(
        ["u_1205", "v_1206", "w_1204", "Werr_1201"], ["E", "N", "U1", "U2"]
    ):
        ds[var] = xr.DataArray(ds["vel"].sel(velbeam=beam), dims=("time", "bindist"))
        ds[var].attrs["units"] = "m s-1"

    ds["u_1205"].attrs["standard_name"] = "eastward_sea_water_velocity"
    ds["v_1206"].attrs["standard_name"] = "northward_sea_water_velocity"
    ds["w_1204"].attrs["standard_name"] = "upward_sea_water_velocity"
    ds["Werr_1201"].attrs[
        "standard_name"
    ] = "indicative_error_from_multibeam_acoustic_doppler_velocity_profiler_in_sea_water"

    ds = aqdutils.magvar_correct(ds)

    ds.attrs["beam_angle"] = ds.attrs["angle"]
    ds.attrs["bin_size"] = ds.attrs["cell"]

    ds = aqdutils.trim_vel(
        ds,
        data_vars=[
            "u_1205",
            "v_1206",
            "w_1204",
            "AGC1_1221",
            "AGC2_1222",
            "AGC3_1223",
            "AGC4_1224",
        ],
    )

    ds = aqdutils.make_bin_depth(ds)

    ds = aqdutils.ds_swap_dims(ds)

    ds = aqdutils.ds_rename(ds)

    ds = qaqc.drop_vars(ds)

    ds = aqdutils.ds_add_attrs(ds, inst_type="RDI")

    ds = ds.drop_vars(["vel", "corr", "int", "Orient", "pg", "velbeam", "beam"])

    nc_filename = Path(ds.attrs["filename"] + "-a.nc")

    if "prefix" in ds.attrs:
        nc_filename = ds.attrs["prefix"] / nc_filename

    ds = utils.check_time_encoding(ds)

    ds.to_netcdf(nc_filename)
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])

    print("Done writing netCDF file", nc_filename)

    return ds
