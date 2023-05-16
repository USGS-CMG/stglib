#!/usr/bin/env python

import yaml

import stglib

args = stglib.cmd.vechdr2cdf_parser().parse_args()

# initialize metadata from the globalatts file
metadata = stglib.read_globalatts(args.gatts)

# Add additional metadata from metadata config file
with open(args.config) as f:
    config = yaml.safe_load(f)

for k in config:
    metadata[k] = config[k]

RAW = stglib.vec.dat2cdf.dat_to_cdf(metadata)
