#!/usr/bin/env python

import stglib

args = stglib.cmd.aqdcdf2nc_parser().parse_args()

if args.atmpres:
    ds = stglib.aqd.cdf2nc.cdf_to_nc(args.cdfname, atmpres=args.atmpres)
else:
    ds = stglib.aqd.cdf2nc.cdf_to_nc(args.cdfname)
