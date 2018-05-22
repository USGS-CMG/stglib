from __future__ import division, print_function
import xarray as xr
from ..core import utils


def cdf_to_nc(cdf_filename, atmpres=None):
    """
    Load raw .cdf file, trim, apply QAQC, and save to .nc
    """

    # Load raw .cdf data
    ds = xr.open_dataset(cdf_filename, autoclose=True)

    # Clip data to in/out water times or via good_ens
    ds = utils.clip_ds(ds)

    if atmpres is not None:
        print("Atmospherically correcting data")

        met = xr.open_dataset(atmpres, autoclose=True)
        # need to save attrs before the subtraction, otherwise they are lost
        # ds['P_1ac'] = ds['P_1'].copy(deep=True)
        attrs = ds['P_1'].attrs
        ds['P_1ac'] = ds['P_1'] - met['atmpres'] - met['atmpres'].offset
        print('Correcting using offset of %f' % met['atmpres'].offset)
        ds['P_1ac'].attrs = attrs

    ds = utils.shift_time(ds,
                          ds.attrs['burst_interval'] *
                          ds.attrs['sample_interval'] / 2)

    ds = utils.create_epic_time(ds)

    ds = ds_add_attrs(ds)

    # add lat/lon coordinates to each variable
    for var in ds.data_vars:
        if 'time' not in var:
            ds = utils.add_lat_lon(ds, var)

    ds = utils.add_min_max(ds)

    ds = utils.add_start_stop_time(ds)

    ds = utils.add_epic_history(ds)

    # Write to .nc file
    print("Writing cleaned/trimmed data to .nc file")
    nc_filename = ds.attrs['filename'] + 'b-cal.nc'

    ds = utils.rename_time(ds)

    ds.to_netcdf(nc_filename, unlimited_dims=['time'])
    print('Done writing netCDF file', nc_filename)

    # rename time variables after the fact to conform with EPIC/CMG standards
    # utils.rename_time(nc_filename)

    return ds

# def da_reshape(ds, var):
#     """
#     Add lon and lat dimensions to DataArrays and reshape to conform to our
#     standard order
#     """
#
#     # Add the dimensions using concat
#     ds[var] = xr.concat([ds[var]], dim=ds['lon'])
#     ds[var] = xr.concat([ds[var]], dim=ds['lat'])
#
#     # Reshape using transpose depending on shape
#     if len(ds[var].shape) == 4:
#         ds[var] = ds[var].transpose('time', 'lon', 'lat', 'frequency')
#     elif len(ds[var].shape) == 3:
#         ds[var] = ds[var].transpose('time', 'lon', 'lat')
#
#     return ds


def ds_add_attrs(ds):
    # Update attributes for EPIC and STG compliance
    ds = utils.ds_coord_no_fillvalue(ds)

    ds['time'].attrs.update({'standard_name': 'time',
                             'axis': 'T'})

    ds['epic_time'].attrs.update({'units': 'True Julian Day',
                                  'type': 'EVEN',
                                  'epic_code': 624})

    ds['epic_time2'].attrs.update({'units': 'msec since 0:00 GMT',
                                   'type': 'EVEN',
                                   'epic_code': 624})

    if 'P_1ac' in ds:
        ds['P_1ac'].attrs.update({'units': 'dbar',
                                  'name': 'Pac',
                                  'long_name': 'Corrected pressure',
                                  '_FillValue': 1e35})
        if 'P_1ac_note' in ds.attrs:
            ds['P_1ac'].attrs.update({'note': ds.attrs['P_1ac_note']})

    return ds
