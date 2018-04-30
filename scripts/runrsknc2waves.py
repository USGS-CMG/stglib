#!/usr/bin/env python

import stglib

args = stglib.cmd.rsknc2waves_parser().parse_args()

ds = stglib.rsk.nc2waves.nc_to_waves(args.ncname)
