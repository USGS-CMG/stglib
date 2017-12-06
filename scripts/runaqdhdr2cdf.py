#!/usr/bin/env python

import sys
sys.path.insert(0, '/Users/dnowacki/Documents/stglib')
import stglib
import argparse
import yaml

parser = argparse.ArgumentParser(description='Convert Aquadopp text files to raw .cdf format. Run this script from the directory containing Aquadopp files')
parser.add_argument('gatts', help='path to global attributes file (gatts formatted)')
parser.add_argument('config', help='path to ancillary config file (YAML formatted)')

args = parser.parse_args()

# initialize metadata from the globalatts file
metadata = stglib.read_globalatts(args.gatts)

# Add additional metadata from metadata config file
config = yaml.safe_load(open(args.config))

for k in config:
    metadata[k] = config[k]

RAW = stglib.aqdhdr2cdf.prf_to_cdf(metadata)
