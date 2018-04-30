#!/usr/bin/env python

import sys
sys.path.insert(0, '/Users/dnowacki/Documents/stglib')
import stglib

args = stglib.cmd.rsknc2waves_parser().parse_args()

ds = stglib.rsk.nc2waves.nc_to_waves(args.ncname)
