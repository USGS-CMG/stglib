#!/usr/bin/env python

import stglib

args = stglib.cmd.sigdolfyncdf2nc_parser().parse_args()

for f in args.cdfname:
    if args.atmpres:
        ds = stglib.sig.dolfyncdf2nc.cdf_to_nc(f, atmpres=args.atmpres)
    else:
        ds = stglib.sig.dolfyncdf2nc.cdf_to_nc(f)
