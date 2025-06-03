import numpy as np
import xarray as xr

from ..core import qaqc, utils
from . import sonutils


def cdf_to_nc(cdf_filename, height=None):
    """
    Load a raw .cdf file and generate a processed .nc file
    """

    # Load raw .cdf data
    ds = xr.open_dataset(cdf_filename)

    # remove units in case we change and we can use larger time steps
    ds.time.encoding.pop("units")

    # Load .nc height above bed data
    ds_hght = xr.open_dataset(height)

    # Reindex bed distance to same timestamp as sonar
    ds_hght = ds_hght.reindex({"time": ds["time"]}, method="nearest")

    instr_bed_dist = ds_hght[ds.attrs["height_var"]]

    # Interpolate NaNs
    instr_bed_dist = instr_bed_dist.interpolate_na(dim="time")

    instr_init_height = ds_hght.attrs["initial_instrument_height"]

    sonar_init_height = ds.attrs["initial_instrument_height"]

    # Correct bed distance for differences in initial height
    hght_diff = instr_init_height - sonar_init_height

    # This works regardless of which instrument is higher because minus a negative is addition
    ds["sonar_hgt"] = instr_bed_dist - hght_diff

    # Add note to sonar height
    ds["sonar_hgt"].attrs.update(
        {
            "note": f"sonar height calculated from {ds.attrs['height_var']} variable in {ds.attrs['brange_file']} using initial height offset of {hght_diff:.2f}"
        }
    )

    # Rename variables (needs to happen before magnetic and theta correction so using new heading name)
    ds = ds_rename_vars(ds)

    # Calculate slant range for each point
    npoints = ds.attrs["SONNDataPoints"]
    total_range = ds.attrs["SONRange"]

    first = total_range / npoints
    # Need to add one more bc of zero indexing
    last = total_range + (total_range / npoints)
    step = total_range / npoints

    # Make slant range for each time point
    ds["SlantRange"] = (
        ("time", "points"),
        np.tile(np.arange(first, last, step), (len(ds.time), 1)),
    )

    # Fill slant range where less than height above bed with NaN
    # > returns where the condition is true and masks everything else
    ds["SlantRange"] = ds["SlantRange"].where(ds["SlantRange"] > ds["sonar_hgt"])

    # Calculate horizontal range (which is rho value needed to plot image)
    ds["HorizontalRange"] = np.sqrt((ds["SlantRange"] ** 2) - (ds["sonar_hgt"] ** 2))

    # Magnetic variation correction (needs to happen before doing theta heading offset corrections)
    ds = magvar_correct(ds)

    # Correct theta and add to ds
    ds["theta"] = correct_theta(ds)

    # Add attributes
    ds = sonutils.ds_add_attrs(ds)

    # QAQC
    ds = qaqc.call_qaqc(ds)

    # Run utilities
    ds = utils.clip_ds(ds)
    ds = utils.create_nominal_instrument_depth(ds)
    ds = utils.create_z(ds)
    ds = utils.ds_add_lat_lon(ds)
    ds = utils.add_start_stop_time(ds)
    ds = utils.add_min_max(ds)
    ds = utils.add_delta_t(ds)
    ds = utils.ds_coord_no_fillvalue(ds)

    # Write to .nc file
    print("Writing cleaned/trimmed data to .nc file")
    nc_filename = f"{ds.attrs['filename']}b_{str(ds.attrs['SONRange'])}m.nc"

    ds.to_netcdf(
        nc_filename, unlimited_dims=["time"], encoding={"time": {"dtype": "i4"}}
    )
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])

    print(f"Done writing netCDF file {nc_filename}")
    #################################################################
    # Average sweeps
    ds_avg = ds.mean(dim="sweep", keep_attrs=True)

    # Recalculate head angle and sonar angle from averaged head position and sonar position
    ds_avg["HeadAngle"].values = 0.3 * (ds_avg["HeadPosition"] - 600)
    ds_avg["SonarAngle"].values = 0.3 * (ds_avg["SonarPosition"] - 600)

    # Vector average headings/pitch/roll over sweeps
    ds_avg["Hdg_1215"].values = vector_avg_angles(ds["Hdg_1215"].values)
    ds_avg["GyroHeading"].values = vector_avg_angles(ds["GyroHeading"].values)
    ds_avg["Ptch_1216"].values = vector_avg_angles(ds["Ptch_1216"].values)
    ds_avg["Roll_1217"].values = vector_avg_angles(ds["Roll_1217"].values)

    # Correct theta
    ds_avg["theta"] = correct_theta(ds_avg)

    # Recalculate min/max
    ds_avg = utils.add_min_max(ds_avg)

    print("Writing averaged data to .nc file")
    nc_avg_filename = f"{ds_avg.attrs['filename']}b_{str(ds.attrs['SONRange'])}m-a.nc"

    ds_avg.to_netcdf(
        nc_avg_filename, unlimited_dims=["time"], encoding={"time": {"dtype": "i4"}}
    )

    utils.check_compliance(nc_avg_filename, conventions=ds_avg.attrs["Conventions"])

    print(f"Finished writing data to {nc_avg_filename}")


def magvar_correct(ds):
    if "Hdg_1215" in ds:

        # Grab out values because can't operate on a variable with more than 2 coordinates
        heading = ds.Hdg_1215.values

        new_heading = heading + ds.attrs["magnetic_variation"].astype(float)
        new_heading = new_heading.round(1)

        # Get heading between 0 and 360
        new_heading = new_heading % 360

        ds["Hdg_1215"].values = new_heading

    return ds


def correct_theta(ds):
    # Reverse image and add 90 degrees to rotate from math to compass convention (north up)
    theta = ds["HeadAngle"]
    theta = -theta + 90

    # Heading correction - Add heading offset to get pointing north
    heading = ds["Hdg_1215"].values
    head_offset = 360 - heading
    theta = theta + head_offset

    # Get theta between 0 and 360
    theta = theta % 360

    # Convert degrees to radians
    theta = np.deg2rad(theta)

    return theta


def ds_rename_vars(ds):
    """
    Rename variables to be CF compliant
    """
    varnames = {"Pitch": "Ptch_1216", "Roll": "Roll_1217", "Heading": "Hdg_1215"}

    # Check to make sure they exist before trying to rename
    for k in varnames:
        if k in ds:
            ds = ds.rename({k: varnames[k]})
    return ds


def vector_avg_angles(angles):

    # Convert angles from degrees to radians
    angles_rad = np.radians(angles)

    # Convert angles to unit vectors
    x = np.cos(angles_rad)
    y = np.sin(angles_rad)

    # Compute the average vector over sweep axis
    avg_x = np.mean(x, axis=1)
    avg_y = np.mean(y, axis=1)

    # Compute the angle of the average vector
    avg_rad = np.arctan2(avg_y, avg_x)

    # Convert the average angle back to degrees
    avg_deg = np.degrees(avg_rad)

    # Get between 0-360
    avg_deg = avg_deg % 360

    return avg_deg
