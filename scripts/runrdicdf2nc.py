#!/usr/bin/env python

import stglib

args = stglib.cmd.rdicdf2nc_parser().parse_args()

if args.atmpres:
    ds = stglib.rdi.cdf2nc.cdf_to_nc(args.cdfname, atmpres=args.atmpres)
else:
    ds = stglib.rdi.cdf2nc.cdf_to_nc(args.cdfname)
