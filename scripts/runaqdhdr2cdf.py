#!/usr/bin/env python

import yaml
from .stglib import eofe
import stglib

args = stglib.cmd.aqdhdr2cdf_parser().parse_args()

# initialize metadata from the globalatts file
metadata = stglib.read_globalatts(args.gatts)

# Add additional metadata from metadata config file
with open(args.config) as f:
    config = yaml.safe_load(f)

for k in config:
    metadata[k] = config[k]

RAW = stglib.aqd.hdr2cdf.prf_to_cdf(metadata)
