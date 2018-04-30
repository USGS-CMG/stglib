#!/usr/bin/env python

import stglib

args = stglib.cmd.exocdf2nc_parser().parse_args()

if args.atmpres:
    ds = stglib.exo.cdf_to_nc(args.cdfname, atmpres=args.atmpres)
else:
    ds = stglib.exo.cdf_to_nc(args.cdfname)
