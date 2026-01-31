import numpy as np
import xarray as xr
from scipy.interpolate import RegularGridInterpolator as RGI
from tqdm import tqdm

from ..core import attrs, qaqc, utils


def nc_to_xy(nc_filename):
    """
    Interpolate sonar data from polar to cartesian coordinates using averaged sweep .nc file
    """

    # Load raw .cdf data
    ds = xr.open_dataset(nc_filename)

    # Convert to x,y coordinates
    final_images = convert_to_xy(ds)

    # Drop variables
    ds = xy_drop_vars(ds)
    ds = ds.drop_dims({"points", "scan"})

    # Add x, y, image to ds
    ds["sonar_image"] = final_images

    # Add attributes
    ds = attrs.ds_add_attrs(ds)

    # QAQC
    ds = qaqc.call_qaqc(ds)

    # Run utilities
    ds = utils.add_min_max(ds)
    ds = utils.ds_coord_no_fillvalue(ds)

    # Write to .nc file
    print("Writing xy data to .nc file")
    nc_filename = f"{ds.attrs['filename']}b_{str(ds.attrs['SONRange'])}m-xy.nc"

    ds.to_netcdf(
        nc_filename, unlimited_dims=["time"], encoding={"time": {"dtype": "i4"}}
    )
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])

    print(f"Done writing netCDF file {nc_filename}")


def xy_drop_vars(ds):
    varnames = list(ds.keys())
    varnames.remove("sonar_hgt")

    for v in varnames:
        if v in ds:
            ds = ds.drop_vars(v)

    return ds


def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)  # rho is the radial distance
    theta = np.arctan2(y, x)  # theta is the angle in radians

    return theta, rho


def convert_to_xy(ds):

    image_list = []
    for j in tqdm(
        np.arange(0, len(ds.time)), bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt}"
    ):

        # Grab variables
        total_range = ds.attrs["SONRange"]
        horz_rng = ds["HorizontalRange"].isel(time=j).values
        theta = ds["theta"].isel(time=j).values
        image = ds["sonar_image"].isel(time=j).values

        # Create x,y grid
        dxy = ds.attrs["dxy"]  # Distance between pixels in resulting image

        x = np.arange(-total_range, total_range + dxy, dxy)
        y = x

        x_grid, y_grid = np.meshgrid(x, y)

        # Convert cartesian grid to polar grid
        theta_grid, rho_grid = cart2pol(x_grid, y_grid)  # radians
        theta_grid = np.degrees(theta_grid)  # Convert radians to degrees
        theta_grid = theta_grid % 360

        # Rotate from math convention to compass convention and get between 0-360 again
        # So rotate 90 degrees from x=0 degrees to y=0 degrees (north up)
        theta_grid = (-theta_grid + 90) % 360

        # Reorder theta and image in ascending order
        idx = np.argsort(theta)
        theta = theta[idx]
        image = image[idx, :]

        # Cut out NaNs in horizontal range because can't interpolate NaNs
        nan_idx = np.where(np.isnan(horz_rng))
        horz_rng = np.delete(horz_rng, nan_idx)
        image = np.delete(image, nan_idx, axis=1)

        theta_deg = np.degrees(theta)

        # Interpolate image with Regular Grid Interpolator
        # Returns a function you have to evaluate
        func = RGI((theta_deg, horz_rng), image, bounds_error=False)

        # Needs extra parentheses because tuple input required
        new_image = func((theta_grid, rho_grid))

        image_list.append(xr.DataArray(new_image, dims=["x", "y"], name="sonar_image"))

    # Concatenate images and add dimensions
    final_images = xr.concat(image_list, dim="time")
    final_images = final_images.assign_coords(x=("x", x))
    final_images = final_images.assign_coords(y=("y", y))

    return final_images
