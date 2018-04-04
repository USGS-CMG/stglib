#!/usr/bin/env python

import sys
sys.path.insert(0, '/Users/dnowacki/Documents/stglib')
import stglib
import argparse

args = stglib.cmd.rskcdf2nc_parse_args()

if args.atmpres:
    ds = stglib.rsk.cdf2nc.cdf_to_nc(args.cdfname, atmpres=args.atmpres)
else:
    ds = stglib.rsk.cdf2nc.cdf_to_nc(args.cdfname)
