import numpy as np
import xarray as xr

from ..core import qaqc, utils


def cdf_to_nc(cdf_filename, atmpres=None, writefile=True, format="NETCDF4"):
    """
    Load raw .cdf file, trim, apply QAQC, and save to .nc
    """

    # Load raw .cdf data
    ds = open_raw_cdf(cdf_filename)

    # Clip data to in/out water times or via good_ens
    ds = utils.clip_ds(ds)

    ds = utils.create_nominal_instrument_depth(ds)

    if atmpres is not None:
        print("Atmospherically correcting data")

        met = xr.load_dataset(atmpres)
        # need to save attrs before the subtraction, otherwise they are lost
        # ds['P_1ac'] = ds['P_1'].copy(deep=True)
        attrs = ds["P_1"].attrs
        ds["P_1ac"] = (
            ds["P_1"]
            - met["atmpres"].reindex_like(
                ds["P_1"], method="nearest", tolerance="10min"
            )
            - met["atmpres"].offset
        )
        print("Correcting using offset of %f" % met["atmpres"].offset)
        ds["P_1ac"].attrs = attrs

    # ds = utils.shift_time(ds,
    #                       ds.attrs['burst_interval'] *
    #                       ds.attrs['sample_interval'] / 2)

    if utils.is_cf(ds):
        pass
    else:
        print("about to create epic times")
        ds = utils.create_epic_times(ds)

        ds = utils.create_2d_time(ds)

    ds = ds_add_attrs(ds)

    if "P_1" in ds:
        ds = ds_add_depth_dim(ds)

    # add lat/lon coordinates to each variable
    # no longer want to do this according to the canonical forms on stellwagen
    # for var in ds.data_vars:
    #     if 'time' not in var:
    #         ds = utils.add_lat_lon(ds, var)

    # trim by minimum pressure for instruments that go out of water_depth
    for v in ["P_1", "P_1ac"]:
        ds = trim_min(ds, v)
        ds = qaqc.trim_bad_ens(ds, v)

    if "Turb" in ds:
        ds = qaqc.trim_min(ds, "Turb")
        ds = qaqc.trim_max(ds, "Turb")
        ds = qaqc.trim_min_diff(ds, "Turb")
        ds = qaqc.trim_max_diff(ds, "Turb")
        ds = qaqc.trim_bad_ens(ds, "Turb")

    ds = utils.add_min_max(ds)

    ds = utils.add_start_stop_time(ds)

    ds = utils.ds_coord_no_fillvalue(ds)

    ds = utils.add_history(ds)

    ds = dw_add_delta_t(ds)

    for var in ds.variables:
        if (var not in ds.coords) and ("time" not in var):
            # cast as float32
            ds = utils.set_var_dtype(ds, var)

    # if we are dealing with continuous instruments, drop sample since it is a singleton dimension
    if "sample" in ds:
        if len(ds["sample"]) == 1:
            ds = ds.squeeze(dim="sample")

    if writefile:
        # Write to .nc file
        print("Writing cleaned/trimmed data to .nc file")
        if "burst" in ds or "sample" in ds:
            nc_filename = ds.attrs["filename"] + "b-cal.nc"
        else:
            nc_filename = ds.attrs["filename"] + "-a.nc"

        ds.to_netcdf(nc_filename, format=format, unlimited_dims=["time"])
        utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])
        print("Done writing netCDF file", nc_filename)

    return ds


def open_raw_cdf(cdf_filename):
    ds = xr.load_dataset(cdf_filename)
    ds.time.encoding.pop(
        "units"
    )  # remove units in case we change and we can use larger time steps
    return ds


def trim_min(ds, var):
    if var + "_min" in ds.attrs:
        print("%s: Trimming using minimum value of %f" % (var, ds.attrs[var + "_min"]))
        # remove full burst if any of the burst values are less than
        # the indicated value
        bads = (ds[var] < ds.attrs[var + "_min"]).any(dim="sample")
        ds[var][bads, :] = np.nan

        notetxt = "Values filled where less than %f units. " % ds.attrs[var + "_min"]

        ds = utils.insert_note(ds, var, notetxt)

    return ds


def ds_add_depth_dim(ds):
    print("Creating depth dimension")
    if "P_1ac" in ds:
        p = "P_1ac"
    else:
        p = "P_1"

    if "NAVD88_ref" in ds.attrs:
        ds["depth"] = xr.DataArray(
            [-ds.attrs["NAVD88_ref"] - ds.attrs["initial_instrument_height"]],
            dims="depth",
        )
        ds["depth"].attrs["VERT_DATUM"] = "NAVD88"
        ds["depth"].attrs["NOTE"] = (
            "Computed as platform depth "
            "[m NAVD88] minus "
            "initial_instrument_height"
        )
    else:
        dim = ["time"]
        if "sample" in ds:
            dim.append("sample")
        ds["depth"] = xr.DataArray(np.atleast_1d(ds[p].mean(dim=dim)), dims="depth")
        ds["depth"].attrs["NOTE"] = "Computed as mean of the pressure sensor"
    ds["depth"].attrs["positive"] = "down"
    ds["depth"].attrs["axis"] = "Z"
    ds["depth"].attrs["units"] = "m"
    ds["depth"].attrs["epic_code"] = 3
    ds["depth"].attrs["standard_name"] = "depth"

    return ds


def ds_add_attrs(ds):
    # Update attributes for EPIC and STG compliance
    ds = utils.ds_coord_no_fillvalue(ds)

    ds["time"].attrs.update(
        {"standard_name": "time", "axis": "T", "long_name": "time (UTC)"}
    )
    ds["time"].encoding["dtype"] = "i4"

    if "sample" in ds:
        ds["sample"].encoding["dtype"] = "i4"
        ds["sample"].attrs["long_name"] = "sample number"
        ds["sample"].attrs["units"] = "1"

    if "epic_time" in ds:
        ds["epic_time"].attrs.update(
            {"units": "True Julian Day", "type": "EVEN", "epic_code": 624}
        )

    if "epic_time2" in ds:
        ds["epic_time2"].attrs.update(
            {"units": "msec since 0:00 GMT", "type": "EVEN", "epic_code": 624}
        )

    if "epic_time_2d" in ds:
        ds["epic_time_2d"].attrs = ds["epic_time"].attrs
    if "epic_time2_2d" in ds:
        ds["epic_time2_2d"].attrs = ds["epic_time2"].attrs

    if "P_1" in ds:
        ds["P_1"].attrs["standard_name"] = "sea_water_pressure"
        ds["P_1"].attrs["long_name"] = "Uncorrected pressure"

    if "P_1ac" in ds:
        ds["P_1ac"].attrs.update(
            {
                "name": "Pac",
                "long_name": "Pressure corrected for changes in atmospheric pressure",
                "standard_name": "sea_water_pressure_due_to_sea_water",
            }
        )
        # ds["P_1ac"].encoding["_FillValue"] = 1e35
        if "P_1ac_note" in ds.attrs:
            ds["P_1ac"].attrs.update({"note": ds.attrs["P_1ac_note"]})

    if "burst" in ds:
        ds["burst"].attrs["units"] = "1"
        # ds["burst"].encoding["_FillValue"] = 1e35
        ds["burst"].attrs["long_name"] = "Burst number"

    if "Turb" in ds:
        ds["Turb"].attrs.update(
            {"long_name": "Turbidity (NTU)", "standard_name": "sea_water_turbidity"}
        )
        # ds["Turb"].encoding["_FillValue"] = 1e35

    if "T_28" in ds:
        ds["T_28"].attrs.update({"standard_name": "sea_water_temperature"})

    if "S_41" in ds:
        ds["S_41"].attrs.update({"standard_name": "sea_water_salinity"})

    if "C_51" in ds:
        ds["C_51"].attrs.update({"standard_name": "sea_water_electrical_conductivity"})

    if "SpC_48" in ds:
        ds["SpC_48"].attrs.update(
            {
                "standard_name": "sea_water_electrical_conductivity",
                "comment": "Temperature compensated to 25 °C",
            }
        )

    if not utils.is_cf(ds):
        ds.attrs["COMPOSITE"] = np.int32(0)
        ds.attrs["COORD_SYSTEM"] = "GEOGRAPHIC + sample"

    return ds


def dw_add_delta_t(ds):

    if "burst_interval" in ds:
        ds.attrs["DELTA_T"] = int(ds.attrs["burst_interval"])

    return ds
