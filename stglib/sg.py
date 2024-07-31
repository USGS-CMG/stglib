import pandas as pd
import xarray as xr

from .core import qaqc, utils


def read_tid(filnam, encoding="utf-8"):
    """Read data from an SBE 26plus Seagauge .tid file into an xarray dataset.

    Parameters
    ----------
    filnam : string
        The filename
    encoding : string, optional
        File encoding. Default 'utf-8'
    Returns
    -------
    xarray.Dataset
        An xarray Dataset of the Seagauge data
    """

    df = pd.read_csv(
        filnam,
        header=None,
        sep=r"\s+",
        names=["Sample", "Date", "Time", "Pressure", "Temp"],
        parse_dates={"time": ["Date", "Time"]},
        encoding=encoding,
        index_col=False,
    )
    df.set_index("time", inplace=True)
    sg = df.to_xarray()
    return sg


def tid_to_cdf(metadata):
    """
    Load a raw .tid file and generate a .cdf file
    """
    basefile = metadata["basefile"]

    ds = read_tid(basefile + ".tid")

    # Convert pressure from psia to dbar
    ds["Pressure"] = ds.Pressure / 14.503773800722 * 10

    ds = utils.write_metadata(ds, metadata)

    ds = utils.ensure_cf(ds)

    cdf_filename = ds.attrs["filename"] + "-raw.cdf"

    ds.to_netcdf(cdf_filename, unlimited_dims=["time"])

    print(f"Finished writing data to {cdf_filename}")

    return ds


def cdf_to_nc(cdf_filename, atmpres=None):
    """
    Load a raw .cdf file and generate a processed .nc file
    """

    # Load raw .cdf data
    ds = xr.open_dataset(cdf_filename)

    # remove units in case we change and we can use larger time steps
    ds.time.encoding.pop("units")

    # Drop sample variable
    ds = ds.drop_vars("Sample")

    # Atmospheric pressure correction
    if atmpres is not None:
        ds = utils.atmos_correct(ds, atmpres)

    # Rename variables to CF compliant names
    ds = ds_rename_vars(ds)

    # Add attributes
    ds = ds_add_attrs(ds)

    # Call QAQC
    ds = sg_qaqc(ds)

    # Run utilities
    ds = utils.create_z(ds)
    ds = utils.clip_ds(ds)
    ds = utils.ds_add_lat_lon(ds)
    ds = utils.create_nominal_instrument_depth(ds)
    ds = utils.add_start_stop_time(ds)
    ds = utils.add_min_max(ds)
    ds = utils.add_delta_t(ds)

    # Write to .nc file
    print("Writing cleaned/trimmed data to .nc file")
    nc_filename = ds.attrs["filename"] + "-a.nc"

    ds.to_netcdf(
        nc_filename, unlimited_dims=["time"], encoding={"time": {"dtype": "i4"}}
    )
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])

    print(f"Done writing netCDF file {nc_filename}")


def ds_rename_vars(ds):
    """
    Rename variables to be CF compliant
    """
    varnames = {"Temp": "T_28", "Pressure": "P_1", "Pressure_ac": "P_1ac"}

    # Check to make sure they exist before trying to rename
    newvars = {}
    for k in varnames:
        if k in ds:
            newvars[k] = varnames[k]
    return ds.rename(newvars)


def ds_add_attrs(ds):
    """
    Add attributes: units, standard name from CF website, long names
    """
    ds = utils.ds_coord_no_fillvalue(ds)

    ds["time"].attrs.update(
        {"standard_name": "time", "axis": "T", "long_name": "time (UTC)"}
    )

    if "T_28" in ds:
        ds["T_28"].attrs.update(
            {
                "units": "degree_C",
                "standard_name": "sea_water_temperature",
                "long_name": "Temperature",
            }
        )

    if "P_1" in ds:
        ds["P_1"].attrs.update(
            {
                "units": "dbar",
                "long_name": "Uncorrected pressure",
                "standard_name": "sea_water_pressure",
            }
        )
    return ds

    if "P_1ac" in ds:
        ds["P_1ac"].attrs.update(
            {
                "units": "dbar",
                "long_name": "Corrected pressure",
                "standard_name": "sea_water_pressure_due_to_seawater",
            }
        )
    return ds


def sg_qaqc(ds):
    """
    QA/QC
    Trim Seagauge data based on metadata
    """

    varlist = ["T_28", "P_1", "P_1ac"]

    [varlist.append(k) for k in ds.data_vars if k not in varlist]

    for var in varlist:
        ds = qaqc.trim_min(ds, var)

        ds = qaqc.trim_max(ds, var)

        ds = qaqc.trim_min_diff(ds, var)

        ds = qaqc.trim_min_diff_pct(ds, var)

        ds = qaqc.trim_max_diff(ds, var)

        ds = qaqc.trim_max_diff_pct(ds, var)

        ds = qaqc.trim_med_diff(ds, var)

        ds = qaqc.trim_med_diff_pct(ds, var)

        ds = qaqc.trim_bad_ens(ds, var)

    for var in varlist:
        ds = qaqc.trim_by_any(
            ds, var
        )  # re-run and trim by other variables as necessary

    return ds


def atmos_correct(ds, atmpres):
    met = xr.load_dataset(atmpres)
    # need to save attrs before the subtraction, otherwise they are lost
    attrs = ds["Pressure"].attrs
    # need to set a tolerance since we can be off by a couple seconds somewhere
    # TODO is this still needed?
    ds["Pressure_ac"] = xr.DataArray(
        ds["Pressure"]
        - met["atmpres"].reindex_like(ds["Pressure"], method="nearest", tolerance="5s")
        - met["atmpres"].attrs["offset"]
    )
    ds["Pressure_ac"].attrs = attrs

    ds.attrs["atmospheric_pressure_correction_file"] = atmpres
    ds.attrs["atmospheric_pressure_correction_offset_applied"] = met["atmpres"].attrs[
        "offset"
    ]

    histtext = f"Atmospherically corrected using time-series from {atmpres} and offset of {met['atmpres'].offset}"

    ds = utils.insert_history(ds, histtext)

    # Also add it as a note
    ds["Pressure_ac"].attrs["note"] = histtext

    return ds
