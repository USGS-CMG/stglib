#!/usr/bin/env python

import stglib

args = stglib.cmd.trollcdf2nc_parser().parse_args()

ds = stglib.troll.cdf_to_nc(args.cdfname)
