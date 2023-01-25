from ..core import utils, waves


def nc_to_waves(nc_filename):

    ds = utils.open_time_2d_dataset(nc_filename)  # this will deal with a cf file, too

    if utils.is_cf(ds):
        pass
    else:
        ds = utils.epic_to_cf_time(ds)

        ds = utils.create_epic_times(ds)

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
    ds = ds_add_attrs(ds)

    # Reshape and associate dimensions with lat/lon
    # if not utils.is_cf:
    #     for var in ["wp_peak", "wh_4061", "wp_4060", "pspec", "water_depth"]:
    #         if var in ds:
    #             ds = utils.add_lat_lon(ds, var)

    # assign min/max (need to do this after trimming):
    ds = utils.add_min_max(ds)

    nc_filename = ds.attrs["filename"] + "s-a.nc"

    ds.to_netcdf(nc_filename, unlimited_dims=["time"])
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])
    print("Done writing netCDF file", nc_filename)

    return ds


def ds_add_attrs(ds):
    """
    Add EPIC and other attributes to variables
    """

    # Update attributes for EPIC and STG compliance
    ds = utils.ds_coord_no_fillvalue(ds)

    ds["time"].attrs.update(
        {"standard_name": "time", "axis": "T", "long_name": "time (UTC)"}
    )

    ds["time"].encoding["dtype"] = "i4"

    if "epic_time" in ds:
        ds["epic_time"].attrs.update(
            {"units": "True Julian Day", "type": "EVEN", "epic_code": 624}
        )

    if "epic_time2" in ds:
        ds["epic_time2"].attrs.update(
            {"units": "msec since 0:00 GMT", "type": "EVEN", "epic_code": 624}
        )

    def add_attributes(var, dsattrs):
        var.attrs.update(
            {
                "initial_instrument_height": dsattrs["initial_instrument_height"],
                "height_depth_units": "m",
            }
        )
        # if "INST_TYPE" in dsattrs:
        #    var.attrs["sensor_type"] = dsattrs["INST_TYPE"]

    ds["wp_peak"].attrs.update(
        {
            "long_name": "Dominant (peak) wave period",
            "units": "s",
            "epic_code": 4063,
            "standard_name": "sea_surface_wave_period_at_variance_spectral_density_maximum",
        }
    )

    ds["wp_4060"].attrs.update(
        {
            "long_name": "Average wave period",
            "units": "s",
            "epic_code": 4060,
            "standard_name": "sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment",
        }
    )

    ds["wh_4061"].attrs.update(
        {
            "long_name": "Significant wave height",
            "units": "m",
            "epic_code": 4061,
            "standard_name": "sea_surface_wave_significant_height",
        }
    )

    ds["pspec"].attrs.update(
        {
            "long_name": "Pressure derived non-directional wave energy spectrum",
            "units": "m^2/Hz",
            "note": "Use caution: all spectra are provisional",
            "standard_name": "sea_surface_wave_variance_spectral_density",
        }
    )

    ds["frequency"].attrs.update({"long_name": "Frequency", "units": "Hz"})

    if "direction" in ds.coords:
        ds["direction"].attrs.update(
            {
                "long_name": "Direction (from, relative to true north)",
                "units": "degrees",
            }
        )

    if "dspec" in ds.data_vars:
        ds["dspec"].attrs.update(
            {
                "long_name": "Directional wave energy spectrum",
                "units": "m^2/Hz/degree",
                "note": "Use caution: all spectra are provisional",
            }
        )

    if "wvdir" in ds.data_vars:
        ds["wvdir"].attrs.update(
            {
                "long_name": (
                    "Direction of peak period " "(from, relative to true north)"
                ),
                "units": "degrees",
                "note": (
                    "Compass direction from which waves are propagating as "
                    "defined by the direction with the greatest energy at "
                    "the peak period"
                ),
            }
        )

    if "dwvdir" in ds.data_vars:
        ds["dwvdir"].attrs.update(
            {
                "long_name": (
                    "Dominant wave direction " "(from, relative to true north)"
                ),
                "units": "degrees",
                "note": (
                    "Compass direction from which waves are propagating as "
                    "defined by the direction band with greatest total "
                    "energy summed over all frequencies"
                ),
            }
        )

    if "wd_4062" in ds.data_vars:
        ds["wd_4062"].attrs.update(
            {
                "long_name": "Mean wave direction",
                "units": "degrees",
                "epic_code": 4062,
                "note": "Compass direction from which waves are propagating",
            }
        )

    for var in [
        "wp_peak",
        "wh_4061",
        "wp_4060",
        "wd_4062",
        "pspec",
        "water_depth",
        "dspec",
        "wvdir",
        "dwvdir",
    ]:
        if var in ds.variables:
            add_attributes(ds[var], ds.attrs)
            ds[var].attrs.update(
                {"minimum": ds[var].min().values, "maximum": ds[var].max().values}
            )

    return ds
