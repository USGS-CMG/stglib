#!/usr/bin/env python

from __future__ import division, print_function

import warnings
import sys
import argparse
import yaml
import xarray as xr
import numpy as np
sys.path.insert(0, '/Users/dnowacki/Documents/rsklib')
sys.path.insert(0, '/Users/dnowacki/Documents/aqdlib')
import rsklib
import aqdlib



def cdf_to_nc(metadata, atmpres=None):
    """
    Load raw .cdf file, trim, apply QAQC, and save to .nc
    """

    cdf_filename = metadata['filename'] + '-raw.cdf'

    ds = xr.open_dataset(cdf_filename, autoclose=True)

    # trim data via one of two methods
    ds = aqdlib.clip_ds(ds, metadata)

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
    for k in ['P_1', 'P_1ac']:
        if k in ds:
            ds[k].attrs.update(minimum=ds[k].min().values, maximum=ds[k].max().values)

            # TODO: published dwave data are not in time, lon, lat, sample format...
            # shouldn't they be?
            # reshape and add lon and lat dimensions

    ds = compute_time(ds)

    ds = ds_add_attrs(ds, metadata)

    ds = rsklib.write_metadata(ds, metadata)

    ds = add_final_metadata(ds)

    # Write to .nc file
    print("Writing cleaned/trimmed data to .nc file")
    write_nc(ds, metadata)

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

def create_epic_time(RAW):

    # create Julian date
    RAW['jd'] = RAW['time'].to_dataframe().index.to_julian_date() + 0.5

    RAW['epic_time'] = np.floor(RAW['jd'])
    if np.all(np.mod(RAW['epic_time'], 1) == 0): # make sure they are all integers, and then cast as such
        RAW['epic_time'] = RAW['epic_time'].astype(np.int32)
    else:
        warnings.warn('not all EPIC time values are integers; '\
                      'this will cause problems with time and time2')

    # TODO: Hopefully this is correct... roundoff errors on big numbers...
    RAW['epic_time2'] = np.round((RAW['jd'] - np.floor(RAW['jd']))*86400000).astype(np.int32)

    return RAW

def write_nc(ds, metadata):
    """Write cleaned and trimmed Dataset to .nc file"""

    nc_filename = metadata['filename'] + 'b-cal.nc'

    ds.to_netcdf(nc_filename, engine='netcdf4')

    # rename time variables after the fact to conform with EPIC/CMG standards
    rsklib.rsknc2diwasp.rename_time(nc_filename)


def add_final_metadata(ds):
    """Add start_time and stop_time global attributes"""

    ds.attrs.update({'start_time': ds['time'][0].values.astype(str),
                     'stop_time': ds['time'][-1].values.astype(str)})

    return ds

def ds_add_attrs(ds, metadata):
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
        if 'P_1ac_note' in metadata:
            ds['P_1ac'].attrs.update({'note': metadata['P_1ac_note']})

    return ds


def main():

    parser = argparse.ArgumentParser(description='Convert raw RBR d|wave .cdf format to processed .nc files')
    parser.add_argument('gatts', help='path to global attributes file (gatts formatted)')
    parser.add_argument('config', help='path to ancillary config file (YAML formatted)')
    parser.add_argument('--atmpres', help='path to cdf file containing atmopsheric pressure data')

    args = parser.parse_args()

    # initialize metadata from the globalatts file
    metadata = rsklib.read_globalatts(args.gatts)

    # Add additional metadata from metadata config file
    config = yaml.safe_load(open(args.config))

    for k in config:
        metadata[k] = config[k]

    if args.atmpres:
        ds = rsklib.cdf_to_nc(metadata, atmpres=args.atmpres)
    else:
        ds = rsklib.cdf_to_nc(metadata)

    return ds

if __name__ == '__main__':
    main()
