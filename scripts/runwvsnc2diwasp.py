#!/usr/bin/env python

import sys
sys.path.insert(0, '/Users/dnowacki/Documents/stglib')
import stglib
import argparse

parser = argparse.ArgumentParser(description='Convert processed Aquadopp waves .nc files using DIWASP')
parser.add_argument('ncname', help='processed .nc filename')

args = parser.parse_args()

ds = stglib.aqd.nc2diwasp.nc_to_diwasp(args.ncname)
