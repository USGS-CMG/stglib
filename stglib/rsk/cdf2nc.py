import numpy as np
import xarray as xr

from ..core import filter, qaqc, utils


def cdf_to_nc(
    cdf_filename, atmpres=None, salwtemp=None, writefile=True, format="NETCDF4"
):
    """
    Load raw .cdf file, trim, apply QAQC, and save to .nc
    """

    # Load raw .cdf data
    ds = open_raw_cdf(cdf_filename)

    is_profile = (
        (ds.attrs["sample_mode"].upper() == "CONTINUOUS")
        and ("featureType" in ds.attrs)
        and (ds.attrs["featureType"] == "profile")
    )

    if is_profile:
        ds = profile_clip_ds(ds)
    else:
        # Clip data to in/out water times or via good_ens
        ds = utils.clip_ds(ds)

        ds = utils.create_nominal_instrument_depth(ds)

    if atmpres is not None and is_profile is False:
        ds = utils.atmos_correct(ds, atmpres)
    elif atmpres is not None and is_profile:
        ds = atmos_correct_profile(ds, atmpres)

    # ds = utils.shift_time(ds,
    #                       ds.attrs['burst_interval'] *
    #                       ds.attrs['sample_interval'] / 2)

    ds = ds_add_attrs(ds, is_profile)

    if not is_profile:
        # add z coordinate dim
        ds = utils.create_z(ds)

        # check for filtered water level
        if "filtered_wl" in ds.attrs and ds.attrs["filtered_wl"].lower() == "true":
            ds = utils.create_water_level_var(ds, salwtemp=salwtemp)
            if "water_level" in ds.data_vars:
                ds = utils.create_filtered_water_level_var(ds)
                ds = ds.drop_vars("water_level")

    # if "P_1" in ds:
    #    ds = ds_add_depth_dim(ds)

    # add lat/lon coordinates to each variable
    # no longer want to do this according to the canonical forms on stellwagen
    # for var in ds.data_vars:
    #     if 'time' not in var:
    #         ds = utils.add_lat_lon(ds, var)

    # trim by minimum pressure for instruments that go out of water_depth
    for v in ["P_1", "P_1ac"]:
        # check if any filtering before other qaqc
        ds = filter.apply_butter_filt(ds, v)
        ds = filter.apply_med_filt(ds, v)

        ds = qaqc.trim_min(ds, v)
        ds = qaqc.trim_bad_ens(ds, v)

    for v in ["Turb", "C_51", "S_41", "T_28", "SpC_48"]:
        if v in ds:
            ds = qaqc.trim_min(ds, v)
            ds = qaqc.trim_max(ds, v)
            ds = qaqc.trim_min_diff(ds, v)
            ds = qaqc.trim_min_diff_pct(ds, v)
            ds = qaqc.trim_max_diff(ds, v)
            ds = qaqc.trim_max_diff_pct(ds, v)
            ds = qaqc.trim_med_diff(ds, v)
            ds = qaqc.trim_med_diff_pct(ds, v)
            ds = qaqc.trim_max_blip(ds, v)
            ds = qaqc.trim_max_blip_pct(ds, v)
            ds = qaqc.trim_bad_ens(ds, v)
            ds = qaqc.trim_bad_ens_indiv(ds, v)
            ds = qaqc.trim_fliers(ds, v)

    # after check for masking vars by other vars
    for var in ds.data_vars:
        ds = qaqc.trim_mask(ds, var)
        ds = qaqc.trim_mask_expr(ds, var)

    if not is_profile:
        """
        # add z coordinate dim
        ds = utils.create_z(ds)

        # check for filtered water level
        if "filtered_wl" in ds.attrs and ds.attrs["filtered_wl"].lower() == "true":
            ds = utils.create_water_level_var(ds, salwtemp=salwtemp)
            ds = utils.create_filtered_water_level_var(ds)
            ds = ds.drop_vars("water_level")
        """

        ds = utils.add_min_max(ds)
        ds = utils.ds_add_lat_lon(ds)

    ds = qaqc.drop_vars(ds)  # Needs to happen after create_water_level
    ds = utils.add_start_stop_time(ds)
    ds = utils.ds_coord_no_fillvalue(ds)
    ds = utils.add_history(ds)
    ds = dw_add_delta_t(ds)

    if is_profile:
        # reset obs and row_start after we are done trimming
        # this is because row_start is supposed to begin at zero according to CF
        # the original obs and row_start are based on the indexes from the raw file as downloaded from the instrument
        # this makes it so there are no skips in obs from removed casts

        attrsbak = ds["obs"].attrs
        obs = np.arange(len(ds["obs"]))
        ds = ds.assign_coords(obs=obs)
        ds["obs"].attrs = attrsbak

        # TODO: this code is mostly redundant with the row_start code in csv2cdf.py. They should be calling the same function
        row_start = np.zeros(ds.row_size.shape, dtype=int)
        for p in range(len(row_start)):
            if p > 0:
                row_start[p] = row_start[p - 1] + ds.row_size[p - 1]
        ds["row_start"].values = row_start

    # if we are dealing with continuous instruments, drop sample since it is a singleton dimension
    if "sample" in ds:
        if len(ds["sample"]) == 1:
            ds = ds.squeeze(dim="sample")

    if writefile:
        # Write to .nc file
        print("Writing cleaned/trimmed data to .nc file")
        if (
            "burst" in ds
            or "sample" in ds
            or ds.attrs["sample_mode"].upper() == "CONTINUOUS"
        ):
            nc_filename = ds.attrs["filename"] + "b.nc"

        elif is_profile:
            nc_filename = ds.attrs["filename"] + "prof.nc"

        else:
            nc_filename = ds.attrs["filename"] + "-a.nc"

        if is_profile:
            ds.to_netcdf(nc_filename, format=format, unlimited_dims=["obs"])
        else:
            ds.to_netcdf(nc_filename, format=format, unlimited_dims=["time"])
        utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])
        print("Done writing netCDF file", nc_filename)

        if (
            "split_profiles" in ds.attrs
            and ds.attrs["split_profiles"].lower() == "true"
        ):
            split_profiles = True
        else:
            split_profiles = False

        if is_profile and split_profiles:
            do_split_profiles(ds)

    return ds


def open_raw_cdf(cdf_filename):
    ds = xr.load_dataset(cdf_filename)
    # remove units in case we change and we can use larger time steps
    ds["time"].encoding.pop("units")
    if "obstime" in ds:
        ds["obstime"].encoding.pop("units")
    return ds


def get_slice(ds, profile):
    row_start = ds.row_start.sel(profile=profile).values
    row_size = ds.row_size.sel(profile=profile).values

    return slice(row_start, row_start + row_size - 1)


def atmos_correct_profile(ds, atmpres):
    met = xr.load_dataset(atmpres)

    # need to save attrs before the subtraction, otherwise they are lost
    attrs = ds["P_1"].attrs
    # apply the correction for each profile in turn. Is there a better way to do this?
    ds["P_1ac"] = xr.full_like(ds["P_1"], np.nan)
    for profile in ds.profile:
        ds["P_1ac"].loc[dict(obs=get_slice(ds, profile))] = (
            ds["P_1"].loc[dict(obs=get_slice(ds, profile))]
            - met["atmpres"].sel(time=ds["time"].sel(profile=profile)).values
            - met["atmpres"].offset
        )
    ds["P_1ac"].attrs = attrs

    ds = utils.insert_history(
        ds,
        f"Atmospherically correcting using time-series from {atmpres} and offset of {met['atmpres'].offset}",
    )
    ds.attrs["atmospheric_pressure_correction_file"] = atmpres
    ds.attrs["atmospheric_pressure_correction_offset_applied"] = met["atmpres"].attrs[
        "offset"
    ]
    if "comment" in met["atmpres"].attrs:
        ds.attrs["atmospheric_pressure_correction_comment"] = met["atmpres"].attrs[
            "comment"
        ]

    return ds


def do_split_profiles(ds):
    max_profile_len = len(str(ds.profile.max().values))
    for profile in ds.profile.values:
        dss = ds.sel(obs=get_slice(ds, profile), profile=profile).copy(deep=True)

        for v in dss.data_vars:
            if "obs" not in dss[v].coords:
                dss[v] = dss[v].expand_dims("profile")  # for CF compliance

        for v in dss.data_vars:
            allnan = True
            if "obs" in dss[v].coords:
                if not np.all(dss[v].isnull()):
                    allnan = False

        dss = utils.insert_history(dss, f"Processed to individual profile #{profile}")

        if allnan:
            print(
                f"All NaN values encountered for profile {profile}; not writing this cast to netCDF"
            )
        else:
            nc_filename = f"{dss.attrs['filename']}prof_{str(profile).zfill(max_profile_len)}-cal.nc"

            # the old unlimited_dims of obs sticks around, so need to specify empty
            dss.to_netcdf(nc_filename, unlimited_dims=[])
            print("Done writing netCDF file", nc_filename)


def trim_min(ds, var):
    if var + "_min" in ds.attrs:
        print(
            "{}: Trimming using minimum value of {:f}".format(
                var, ds.attrs[var + "_min"]
            )
        )
        # remove full burst if any of the burst values are less than
        # the indicated value

        if "sample" in ds:
            bads = (ds[var] < ds.attrs[var + "_min"]).any(dim="sample")
            ds[var][bads, :] = np.nan
        else:
            ds[var][ds[var] < ds.attrs[var + "_min"]] = np.nan

        notetxt = "Values filled where less than %f units. " % ds.attrs[var + "_min"]

        ds = utils.insert_note(ds, var, notetxt)

    return ds


# def ds_add_depth_dim(ds):
#     print("Creating depth dimension")
#     if "P_1ac" in ds:
#         p = "P_1ac"
#     else:
#         p = "P_1"
#
#     if "NAVD88_ref" in ds.attrs:
#         ds["depth"] = xr.DataArray(
#             [-ds.attrs["NAVD88_ref"] - ds.attrs["initial_instrument_height"]],
#             dims="depth",
#         )
#         ds["depth"].attrs["VERT_DATUM"] = "NAVD88"
#         ds["depth"].attrs["NOTE"] = (
#             "Computed as platform depth "
#             "[m NAVD88] minus "
#             "initial_instrument_height"
#         )
#     else:
#         dim = ["time"]
#         if "sample" in ds:
#             dim.append("sample")
#         ds["depth"] = xr.DataArray(np.atleast_1d(ds[p].mean(dim=dim)), dims="depth")
#         ds["depth"].attrs["NOTE"] = "Computed as mean of the pressure sensor"
#     ds["depth"].attrs["positive"] = "down"
#     ds["depth"].attrs["axis"] = "Z"
#     ds["depth"].attrs["units"] = "m"
#     ds["depth"].attrs["epic_code"] = 3
#     ds["depth"].attrs["standard_name"] = "depth"
#
#     return ds


def ds_add_attrs(ds, is_profile):
    # Update attributes for EPIC and STG compliance
    ds = utils.ds_coord_no_fillvalue(ds)

    ds["time"].attrs.update(
        {"standard_name": "time", "axis": "T", "long_name": "time (UTC)"}
    )
    if is_profile:
        ds["time"].attrs["long_name"] = "observation time (UTC)"

    if "sample" in ds:
        ds["sample"].attrs["long_name"] = "sample number"
        ds["sample"].attrs["units"] = "1"

    if "P_1" in ds:
        ds["P_1"].attrs["standard_name"] = "sea_water_pressure"
        ds["P_1"].attrs["long_name"] = "Uncorrected pressure"
        ds["P_1"].attrs["epic_code"] = 1

    if "P_1ac" in ds:
        ds["P_1ac"].attrs.update(
            {
                "long_name": "Corrected pressure",
                "standard_name": "sea_water_pressure_due_to_sea_water",
            }
        )
        if "P_1ac_note" in ds.attrs:
            ds["P_1ac"].attrs.update({"note": ds.attrs["P_1ac_note"]})

    if "burst" in ds:
        ds["burst"].attrs["units"] = "1"
        ds["burst"].attrs["long_name"] = "Burst number"

    if "Turb" in ds:
        ds["Turb"].attrs.update(
            {"long_name": "Turbidity (NTU)", "standard_name": "sea_water_turbidity"}
        )

    if "T_28" in ds:
        ds["T_28"].attrs.update({"standard_name": "sea_water_temperature"})

    if "S_41" in ds:
        ds["S_41"].attrs.update({"standard_name": "sea_water_practical_salinity"})

    if "C_51" in ds:
        ds["C_51"].attrs.update({"standard_name": "sea_water_electrical_conductivity"})

    if "SpC_48" in ds:
        ds["SpC_48"].attrs.update(
            {
                "standard_name": "sea_water_electrical_conductivity_at_reference_temperature",
                "comment": "Temperature compensated to 25 °C",
            }
        )

    def add_attributes(var, dsattrs):
        var.attrs.update(
            {
                "initial_instrument_height": dsattrs["initial_instrument_height"],
                "height_depth_units": "m",
            }
        )

    # for var in ds.variables:
    #    if (var not in ds.coords) and ("time" not in var):
    #        add_attributes(ds[var], ds.attrs)

    return ds


def dw_add_delta_t(ds):
    if "burst_interval" in ds.attrs:
        ds.attrs["DELTA_T"] = int(ds.attrs["burst_interval"])

    return ds


def profile_clip_ds(ds):
    print(
        f"first profile in full file: {ds['time'].min().values}, idx {np.argmin(ds['time'].values)}"
    )
    print(
        f"last profile in full file: {ds['time'].max().values}, idx {np.argmax(ds['time'].values)}"
    )
    if "good_ens" in ds.attrs:
        # we have good ensemble indices in the metadata

        # so we can deal with multiple good_ens ranges, or just a single range
        # these are good profile numbers
        good_ens = ds.attrs["good_ens"]
        goods = []
        for n in range(0, len(good_ens), 2):
            goods.append(np.arange(good_ens[n], good_ens[n + 1]))
        goods = np.hstack(goods)

        # these are good obs numbers
        goodobs = []
        for profile in ds.profile.values:
            if profile in goods:
                s = get_slice(ds, profile)
                goodobs.append(list(range(s.start, s.stop + 1)))
        goodobs = np.hstack(goodobs)

        ds = ds.sel(profile=goods, obs=goodobs)

        histtext = "Data clipped using good_ens values of {}.".format(str(good_ens))

        ds = utils.insert_history(ds, histtext)

    else:
        print("Did not clip data; no values specified in metadata")

    print(
        f"first profile in clipped file: {ds['time'].min().values}, idx {np.argmin(ds['time'].values)}"
    )
    print(
        f"last profile in clipped file: {ds['time'].max().values}, idx {np.argmax(ds['time'].values)}"
    )

    return ds
