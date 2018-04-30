#!/usr/bin/env python

import stglib
import argparse

args = stglib.cmd.rskcdf2nc_parser().parse_args()

if args.atmpres:
    ds = stglib.rsk.cdf2nc.cdf_to_nc(args.cdfname, atmpres=args.atmpres)
else:
    ds = stglib.rsk.cdf2nc.cdf_to_nc(args.cdfname)
