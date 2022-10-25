#!/usr/bin/env python

import stglib

args = stglib.cmd.veccdf2nc_parser().parse_args()

if args.atmpres:
    ds = stglib.vec.cdf2nc.cdf_to_nc(args.cdfname, atmpres=args.atmpres)
else:
    ds = stglib.vec.cdf2nc.cdf_to_nc(args.cdfname)
