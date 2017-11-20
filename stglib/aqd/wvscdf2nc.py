#!/usr/bin/env python

from __future__ import division, print_function
import sys
sys.path.insert(0, '/Users/dnowacki/Documents/aqdlib')
import aqdlib
from aqdhdr2cdf import compute_time, read_aqd_hdr, check_metadata, check_orientation
import aqdcdf2nc
import qaqc

def cdf_to_nc(cdf_filename, metadata, atmpres=False):

    # Load raw .cdf data
    VEL = aqdcdf2nc.load_cdf(cdf_filename, metadata, atmpres=atmpres)

    # Clip data to in/out water times or via good_ens
    VEL = aqdcdf2nc.clip_ds(VEL, metadata)

    # Create water_depth variables
    VEL, metadata = qaqc.create_water_depth(VEL, metadata)

    # Create depth variable depending on orientation
    VEL, T = qaqc.set_orientation(VEL, VEL['TransMatrix'].values, metadata)

    VEL = qaqc.make_bin_depth(VEL, metadata)

    VEL = aqdcdf2nc.ds_rename(VEL, waves=True)

    VEL = aqdcdf2nc.ds_add_attrs(VEL, metadata, waves=True)

    # TODO: Need to add all global attributes from CDF to NC file (or similar)
    VEL = qaqc.add_min_max(VEL)

    nc_filename = metadata['filename'] + 'wvsb-cal.nc' # TODO: why is a "b" in there?
    VEL.to_netcdf(nc_filename, unlimited_dims='time')
    print('Done writing netCDF file', nc_filename)

    # rename time variables after the fact to conform with EPIC/CMG standards
    aqdlib.aqdcdf2nc.rename_time(nc_filename)

    return VEL


def main():
    import sys
    sys.path.insert(0, '/Users/dnowacki/Documents/aqdlib')
    import aqdlib
    import argparse
    import yaml

    parser = argparse.ArgumentParser(description='Convert raw Aquadopp .cdf wave files to processed .nc files')
    parser.add_argument('cdfname', help='raw .CDF filename')
    parser.add_argument('gatts', help='path to global attributes file (gatts formatted)')
    parser.add_argument('config', help='path to ancillary config file (YAML formatted)')
    parser.add_argument('--atmpres', help='path to cdf file containing atmopsheric pressure data')

    args = parser.parse_args()

    # initialize metadata from the globalatts file
    metadata = aqdlib.read_globalatts(args.gatts)

    # Add additional metadata from metadata config file
    config = yaml.safe_load(open(args.config))

    for k in config:
        metadata[k] = config[k]

    if args.atmpres:
        VEL = cdf_to_nc(args.cdfname, metadata, atmpres=args.atmpres)
    else:
        VEL = cdf_to_nc(args.cdfname, metadata)

if __name__ == '__main__':
    main()
