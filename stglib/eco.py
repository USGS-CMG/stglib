import numpy as np
import pandas as pd
import xarray as xr

from .core import qaqc, utils


def read_par(filnam, spb=False, skiprows=None, skipfooter=0):
    """Read data from a WET Labs PAR csv file into an xarray
    Dataset.

    Parameters
    ----------
    filnam : string
        The filename
    spb: bool, optional
        Samples per burst if using burst sampling
    skiprows : int, optional
        How many header rows to skip. Default None
    skipfooter : into, optional
        How many footer rows to skip. Default 0

    Returns
    -------
    xarray.Dataset
        An xarray Dataset of the PAR data
    """

    names = ["date", "time", "counts"]

    par = read_eco_csv(filnam, names, skiprows=skiprows, skipfooter=skipfooter)

    return eco_pd_to_xr(par, spb=spb)


def read_ntu(filnam, spb=False, skiprows=None, skipfooter=0):
    """Read data from a WET Labs NTU csv file into an xarray
    Dataset.

    Parameters
    ----------
    filnam : string
        The filename
    spb: bool, optional
        Samples per burst if using burst sampling
    skiprows : int, optional
        How many header rows to skip. Default None
    skipfooter : into, optional
        How many footer rows to skip. Default 0

    Returns
    -------
    xarray.Dataset
        An xarray Dataset of the PAR data
    """

    names = ["date", "time", "a", "counts", "b"]

    ntu = read_eco_csv(filnam, names, skiprows=skiprows, skipfooter=skipfooter)

    return eco_pd_to_xr(ntu, spb=spb)


def read_eco_csv(filnam, names, skiprows=None, skipfooter=0):

    return pd.read_csv(
        filnam,
        sep="\t",
        names=names,
        parse_dates=[["date", "time"]],
        infer_datetime_format=True,
        engine="python",
        skiprows=skiprows,
        skipfooter=skipfooter,
    )


def eco_pd_to_xr(df, spb=False):

    if spb:
        # get middle time
        times = df["date_time"].values.reshape((-1, spb))[:, int(spb / 2)]
        counts = df["counts"].values.reshape((-1, spb))
        sample = range(spb)

        ds = xr.Dataset(
            {
                "time": ("time", times),
                "counts": (["time", "sample"], counts),
                "sample": ("sample", sample),
            }
        )
    else:
        times = df["date_time"]
        counts = df["counts"]

        ds = xr.Dataset({"time": ("time", times), "counts": ("time", counts)})

    return ds


def csv_to_cdf(metadata):
    """
    Process ECO .csv file to a raw .cdf file
    """

    basefile = metadata["basefile"]

    if "par" in metadata["INST_TYPE"].lower():
        f = read_par
    elif "ntu" in metadata["INST_TYPE"].lower():
        f = read_ntu
    kwargs = {
        "spb": metadata["spb"],
        "skiprows": metadata["skiprows"],
        "skipfooter": metadata["skipfooter"],
    }
    try:
        ds = f(basefile, **kwargs)
    except UnicodeDecodeError:
        # try reading as Mac OS Western for old versions of Mac Excel
        ds = f(basefile, encoding="mac-roman", **kwargs)

    metadata.pop("skiprows")
    metadata.pop("skipfooter")

    # write out metadata first, then deal exclusively with xarray attrs
    ds = utils.write_metadata(ds, metadata)

    del metadata

    ds = utils.ensure_cf(ds)

    # configure file
    cdf_filename = ds.attrs["filename"] + "-raw.cdf"

    ds.to_netcdf(cdf_filename, unlimited_dims=["time"])

    print("Finished writing data to %s" % cdf_filename)

    return ds


def cdf_to_nc(cdf_filename, atmpres=False):
    """
    Load a "raw" .cdf file and generate a processed .nc file
    """

    # Load raw .cdf data
    ds = xr.open_dataset(cdf_filename)

    # definition of PAR is
    # PAR = Im * 10 ^ ((x-a0)/a1)
    # Where
    # Im is the immersion coefficient
    # a1 is the scaling factor
    # a0 is the voltage offset, typically 0
    # x is the voltage
    # The manufacturer calculates PAR in units of μmol photons/m2/s1
    # from Sea-Bird Scientific, ECO PAR User Manual
    # Document No. par170706, 2017-07-06, Version B
    # https://www.seabird.com/asset-get.download.jsa?id=54627862518

    if "par" in ds.attrs["INST_TYPE"].lower():
        ds["PAR_905"] = ds.attrs["Im"] * 10 ** (
            (ds["counts"].mean(dim="sample") - ds.attrs["a0"]) / ds.attrs["a1"]
        )
        ds["PAR_905"].attrs["units"] = "umol m-2 s-1"
        ds["PAR_905"].attrs["long_name"] = "Photosynthetically active " "radiation"

    if "ntu" in ds.attrs["INST_TYPE"].lower():
        if "user_ntucal_coeffs" in ds.attrs:
            ds["Turb"] = xr.DataArray(
                np.polyval(ds.attrs["user_ntucal_coeffs"], ds["counts"]),
                dims=["time", "sample"],
            ).mean(dim="sample")
            ds["Turb"].attrs["units"] = "1"
            ds["Turb"].attrs["long_name"] = "Turbidity (NTU)"
            ds["Turb"].attrs["standard_name"] = "sea_water_turbidity"
            ds["Turb_std"] = xr.DataArray(
                np.polyval(ds.attrs["user_ntucal_coeffs"], ds["counts"]),
                dims=["time", "sample"],
            ).std(dim="sample")
            ds["Turb_std"].attrs["units"] = "1"
            ds["Turb_std"].attrs[
                "long_name"
            ] = "Turbidity burst standard deviation (NTU)"
            ds["Turb_std"].attrs["standard_name"] = "sea_water_turbidity"
            ds["Turb_std"].attrs["cell_methods"] = "time: standard_deviation"

    ds = ds.drop(["counts", "sample"])

    # Clip data to in/out water times or via good_ens
    ds = utils.clip_ds(ds)

    ds = eco_qaqc(ds)

    # assign min/max:
    ds = utils.add_min_max(ds)

    ds = utils.add_start_stop_time(ds)

    ds = utils.add_delta_t(ds)

    # add lat/lon coordinates
    ds = utils.ds_add_lat_lon(ds)

    ds = ds_add_attrs(ds)

    ds = utils.create_z(ds)

    # add lat/lon coordinates to each variable
    # for var in ds.variables:
    #     if (var not in ds.coords) and ("time" not in var):
    #         # ds = utils.add_lat_lon(ds, var)
    #         # ds = utils.no_p_add_depth(ds, var)
    #         ds = utils.add_z_if_no_pressure(ds, var)
    #         # cast as float32
    #         # ds = utils.set_var_dtype(ds, var)

    # Write to .nc file
    print("Writing cleaned/trimmed data to .nc file")
    nc_filename = ds.attrs["filename"] + "-a.nc"

    ds.to_netcdf(
        nc_filename, unlimited_dims=["time"], encoding={"time": {"dtype": "i4"}}
    )
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])
    print("Done writing netCDF file", nc_filename)


def ds_add_attrs(ds):
    # Update attributes for EPIC and STG compliance
    ds = utils.ds_coord_no_fillvalue(ds)

    ds["time"].attrs.update(
        {"standard_name": "time", "axis": "T", "long_name": "time (UTC)"}
    )

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
                # 'nominal_instrument_depth': dsattrs['nominal_instrument_depth'],
                "height_depth_units": "m",
            }
        )

    for var in ds.variables:
        if (var not in ds.coords) and ("time" not in var):
            add_attributes(ds[var], ds.attrs)

    return ds


def eco_qaqc(ds):
    # QA/QC ECO data
    if "ntu" in ds.attrs["INST_TYPE"].lower():
        for var in ["Turb"]:
            ds = qaqc.trim_min_diff(ds, var)

            ds = qaqc.trim_max_diff(ds, var)

            ds = qaqc.trim_med_diff(ds, var)

            ds = qaqc.trim_med_diff_pct(ds, var)

            ds = qaqc.trim_maxabs_diff_2d(ds, var)

            ds = qaqc.trim_maxabs_diff(ds, var)

            ds = qaqc.trim_max_std(ds, var)

            ds = qaqc.trim_min(ds, var)

            ds = qaqc.trim_max(ds, var)

            ds = qaqc.trim_bad_ens(ds, var)

        # after check for masking vars by others
        for var in ["Turb"]:
            ds = qaqc.trim_mask(ds, var)

    return ds
