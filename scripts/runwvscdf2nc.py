#!/usr/bin/env python

import stglib

args = stglib.cmd.wvscdf2nc_parser().parse_args()

if args.atmpres:
    VEL = stglib.aqd.wvscdf2nc.cdf_to_nc(args.cdfname, atmpres=args.atmpres)
else:
    VEL = stglib.aqd.wvscdf2nc.cdf_to_nc(args.cdfname)
