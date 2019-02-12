#!/usr/bin/env python

import stglib

args = stglib.cmd.hlmcdf2nc_parser().parse_args()

ds = stglib.hobo.cdf_to_nc(args.cdfname)
