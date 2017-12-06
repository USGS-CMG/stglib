#!/usr/bin/env python

import sys
sys.path.insert(0, '/Users/dnowacki/Documents/stglib')
import stglib
import argparse
import yaml

parser = argparse.ArgumentParser(description='Convert raw Aquadopp .cdf format to processed .nc files')
parser.add_argument('cdfname', help='raw .CDF filename')
parser.add_argument('gatts', help='path to global attributes file (gatts formatted)')
parser.add_argument('config', help='path to ancillary config file (YAML formatted)')
parser.add_argument('--atmpres', help='path to cdf file containing atmopsheric pressure data')

args = parser.parse_args()

# initialize metadata from the globalatts file
metadata = stglib.read_globalatts(args.gatts)

# Add additional metadata from metadata config file
config = yaml.safe_load(open(args.config))

for k in config:
    metadata[k] = config[k]

if args.atmpres:
    VEL = stglib.aqdcdf2nc.cdf_to_nc(args.cdfname, metadata, atmpres=args.atmpres)
else:
    VEL = stglib.aqdcdf2nc.cdf_to_nc(args.cdfname, metadata)
