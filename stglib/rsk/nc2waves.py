from ..core import utils, waves


def nc_to_waves(nc_filename):

    ds = utils.open_time_2d_dataset(nc_filename)  # this will deal with a cf file, too

    spec = waves.make_waves_ds(ds)

    for k in ["wp_peak", "wh_4061", "wp_4060", "pspec"]:
        ds[k] = spec[k]

    # ds = utils.create_water_depth(ds)

    ds = utils.create_water_depth_var(ds)

    for k in ["P_1", "P_1ac", "sample", "T_28"]:
        if k in ds:
            ds = ds.drop_vars(k)

    ds = utils.trim_max_wp(ds)

    ds = utils.trim_min_wh(ds)

    ds = utils.trim_max_wh(ds)

    ds = utils.trim_wp_ratio(ds)

    # Add attrs
    ds = utils.ds_add_wave_attrs(ds)

    # assign min/max (need to do this after trimming):
    ds = utils.add_min_max(ds)

    nc_filename = ds.attrs["filename"] + "s-a.nc"

    ds.to_netcdf(nc_filename, unlimited_dims=["time"])
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])
    print("Done writing netCDF file", nc_filename)

    return ds
