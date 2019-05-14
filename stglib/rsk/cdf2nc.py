from __future__ import division, print_function
import numpy as np
import xarray as xr
from ..core import utils


def cdf_to_nc(cdf_filename,
              atmpres=None,
              writefile=True,
              format='NETCDF3_64BIT'):
    """
    Load raw .cdf file, trim, apply QAQC, and save to .nc
    """

    # Load raw .cdf data
    ds = xr.open_dataset(cdf_filename).load()
    ds.close()

    # Clip data to in/out water times or via good_ens
    ds = utils.clip_ds(ds)

    ds = utils.create_nominal_instrument_depth(ds)

    if atmpres is not None:
        print("Atmospherically correcting data")

        met = xr.open_dataset(atmpres).load()
        met.close()
        # need to save attrs before the subtraction, otherwise they are lost
        # ds['P_1ac'] = ds['P_1'].copy(deep=True)
        attrs = ds['P_1'].attrs
        ds['P_1ac'] = ds['P_1'] - met['atmpres'] - met['atmpres'].offset
        print('Correcting using offset of %f' % met['atmpres'].offset)
        ds['P_1ac'].attrs = attrs

    # ds = utils.shift_time(ds,
    #                       ds.attrs['burst_interval'] *
    #                       ds.attrs['sample_interval'] / 2)

    ds = utils.create_epic_times(ds)

    ds = utils.create_2d_time(ds)

    ds = ds_add_attrs(ds)

    ds = ds_add_depth_dim(ds)

    # add lat/lon coordinates to each variable
    # no longer want to do this according to the canonical forms on stellwagen
    # for var in ds.data_vars:
    #     if 'time' not in var:
    #         ds = utils.add_lat_lon(ds, var)

    # trim by minimum pressure for instruments that go out of water_depth
    ds = trim_min(ds, 'P_1ac')

    ds = utils.add_min_max(ds)

    ds = utils.add_start_stop_time(ds)

    ds = utils.ds_coord_no_fillvalue(ds)

    ds = utils.add_history(ds)

    ds = dw_add_delta_t(ds)

    for var in ds.variables:
        if (var not in ds.coords) and ('time' not in var):
            # cast as float32
            ds = utils.set_var_dtype(ds, var)

    if writefile:
        # Write to .nc file
        print("Writing cleaned/trimmed data to .nc file")
        nc_filename = ds.attrs['filename'] + 'b-cal.nc'

        ds.to_netcdf(nc_filename, format=format, unlimited_dims=['time'])

        # Rename time variables for EPIC compliance, keeping a time_cf
        # coorindate.
        utils.rename_time_2d(nc_filename, ds)

        print('Done writing netCDF file', nc_filename)

    # rename time variables after the fact to conform with EPIC/CMG standards
    # utils.rename_time(nc_filename)

    return ds

def trim_min(ds, var):
    if var + '_min' in ds.attrs:
        print('%s: Trimming using minimum value of %f' %
              (var, ds.attrs[var + '_min']))
        # remove full burst if any of the burst values are less than
        # the indicated value
        bads = (ds[var] < ds.attrs[var + '_min']).any(dim='sample')
        ds[var][bads, :] = np.nan

        notetxt = ('Values filled where less than %f units. ' %
                   ds.attrs[var + '_min'])

        ds = utils.insert_note(ds, var, notetxt)

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

def ds_add_depth_dim(ds):
    print('Creating depth dimension')
    if 'P_1ac' in ds:
        p = 'P_1ac'
    else:
        p = 'P_1'

    if 'NAVD88_ref' in ds.attrs:
        ds['depth'] = xr.DataArray([-ds.attrs['NAVD88_ref'] -
            ds.attrs['initial_instrument_height']], dims='depth')
        ds['depth'].attrs['VERT_DATUM'] = 'NAVD88'
        ds['depth'].attrs['NOTE'] = ('Computed as platform depth '
                                     '[m NAVD88] minus '
                                     'initial_instrument_height')
    else:
        ds['depth'] = xr.DataArray([ds[p].mean(dim=['time', 'sample'])],
                                   dims='depth')
        ds['depth'].attrs['NOTE'] = 'Computed as mean of the pressure sensor'
    ds['depth'].attrs['positive'] = 'down'
    ds['depth'].attrs['axis'] = 'z'
    ds['depth'].attrs['units'] = 'm'
    ds['depth'].attrs['epic_code'] = 3

    return ds


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

    if 'epic_time_2d' in ds:
        ds['epic_time_2d'].attrs = ds['epic_time'].attrs
    if 'epic_time2_2d' in ds:
        ds['epic_time2_2d'].attrs = ds['epic_time2'].attrs

    if 'P_1ac' in ds:
        ds['P_1ac'].attrs.update({'units': 'dbar',
                                  'name': 'Pac',
                                  'long_name': 'Corrected pressure'})
        ds['P_1ac'].encoding['_FillValue'] = 1e35
        if 'P_1ac_note' in ds.attrs:
            ds['P_1ac'].attrs.update({'note': ds.attrs['P_1ac_note']})

    if 'burst' in ds:
        ds['burst'].encoding['_FillValue'] = 1e35

    ds.attrs['COMPOSITE'] = np.int32(0)
    ds.attrs['COORD_SYSTEM'] = 'GEOGRAPHIC + sample'

    return ds

def dw_add_delta_t(ds):

    ds.attrs['DELTA_T'] = int(ds.attrs['burst_interval'])

    return ds
