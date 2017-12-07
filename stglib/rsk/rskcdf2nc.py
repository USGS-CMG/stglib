from __future__ import division, print_function

import warnings
import xarray as xr
import numpy as np
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

    # assign min/max:
    ds = utils.add_min_max(ds)

    ds = compute_time(ds)

    ds = ds_add_attrs(ds)

    ds = add_final_rsk_metadata(ds)

    # Write to .nc file
    print("Writing cleaned/trimmed data to .nc file")
    nc_filename = ds.attrs['filename'] + 'b-cal.nc'

    ds.to_netcdf(nc_filename)
    print('Done writing netCDF file', nc_filename)

    # rename time variables after the fact to conform with EPIC/CMG standards
    utils.rename_time(nc_filename)

    return ds


def compute_time(RAW):
    """Compute Julian date and then time and time2 for use in netCDF file"""

    # shift times to center of ensemble
    timeshift = RAW.attrs['burst_interval']*RAW.attrs['sample_interval']/2

    if timeshift.is_integer():
        RAW['time'] = RAW['time'] + np.timedelta64(int(timeshift), 's')
        print('Time shifted by:', int(timeshift), 's')
    else:
        warnings.warn('time NOT shifted because not a whole number of seconds: %f s ***' % timeshift)

    RAW = create_epic_time(RAW)

    return RAW

def add_final_rsk_metadata(ds):
    """Add start_time and stop_time global attributes"""

    ds.attrs['history'] = 'Processed to EPIC using rskcdf2nc.py. ' + ds.attrs['history']

    ds.attrs.update({'start_time': ds['time'][0].values.astype(str),
                     'stop_time': ds['time'][-1].values.astype(str)})

    return ds

def ds_add_attrs(ds):
    # Update attributes for EPIC and STG compliance
    ds.lat.encoding['_FillValue'] = False
    ds.lon.encoding['_FillValue'] = False
    ds.depth.encoding['_FillValue'] = False
    ds.time.encoding['_FillValue'] = False
    ds.epic_time.encoding['_FillValue'] = False
    ds.epic_time2.encoding['_FillValue'] = False
    ds.sample.encoding['_FillValue'] = False

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
