#!/usr/bin/env python

import sys
sys.path.insert(0, '/Users/dnowacki/Documents/stglib')
import stglib

args = stglib.cmd.rsknc2diwasp_parse_args()

ds = stglib.rsk.nc2diwasp.nc_to_diwasp(args.ncname)
