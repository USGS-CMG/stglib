#!/usr/bin/env python

import stglib

args = stglib.cmd.wvsnc2waves_parser().parse_args()

ds = stglib.aqd.wvsnc2waves.nc_to_waves(args.ncname)
