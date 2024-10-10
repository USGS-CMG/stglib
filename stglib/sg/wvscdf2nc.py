import xarray as xr

from ..core import utils
from . import sgutils


def cdf_to_nc(cdf_filename, atmpres=None):
    """
    Load a raw .cdf file and generate a processed .nc file
    """

    # Check for atmpres correction file
    # Atmpres file is required for seagauge because pressure is measured as absolute pressure
    if atmpres is None:
        raise FileNotFoundError(
            "The atmpres file does not exist. Atmpres file is required for Seagauge because pressure is measured as absolute pressure."
        )

    # Load raw .cdf data
    ds = xr.open_dataset(cdf_filename)

    # remove units in case we change and we can use larger time steps
    ds.time.encoding.pop("units")

    # Add sample_interval to metadata (Convert Hertz to sample interval in seconds)
    ds.attrs["sample_interval"] = 1 / float(ds.attrs["sample_rate"])

    # Atmospheric pressure correction
    ds = sgutils.atmos_correct_burst(ds, atmpres)

    # Drop variables
    ds = ds.drop("burst_number")

    # Edit metadata depending
    ds = ds_drop_meta(ds)

    # Add attributes
    ds = sgutils.ds_add_attrs(ds)

    # Call QAQC
    ds = sgutils.sg_qaqc(ds)

    # Run utilities
    ds = utils.clip_ds(ds)
    ds = utils.create_nominal_instrument_depth(ds)
    ds = utils.create_z(ds)
    ds = utils.ds_add_lat_lon(ds)
    ds = utils.add_start_stop_time(ds)
    ds = utils.add_min_max(ds)

    # Write to .nc file
    print("Writing cleaned/trimmed data to .nc file")
    nc_filename = ds.attrs["filename"] + "b-cal.nc"

    ds.to_netcdf(
        nc_filename, unlimited_dims=["time"], encoding={"time": {"dtype": "i4"}}
    )
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])

    print(f"Done writing netCDF file {nc_filename}")


def ds_drop_meta(ds):
    """
    Drop global attribute metadata not needed for .wb file
    """
    gatts = [
        "TideInterval",
        "TideIntervalUnits",
        "TideDuration",
        "TideDurationUnits",
        "TideSamplesPerDay",
        "NumberOfTideMeasurements",
    ]

    # Check to make sure they exist
    for k in gatts:
        if k in ds.attrs:
            del ds.attrs[k]
    return ds
