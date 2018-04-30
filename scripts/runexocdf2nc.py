#!/usr/bin/env python

import stglib
import argparse

parser = argparse.ArgumentParser(description='Convert raw EXO .cdf format to processed .nc files')
parser.add_argument('cdfname', help='raw .CDF filename')
parser.add_argument('--atmpres', help='path to cdf file containing atmopsheric pressure data')

args = parser.parse_args()

if args.atmpres:
    ds = stglib.exo.cdf_to_nc(args.cdfname, atmpres=args.atmpres)
else:
    ds = stglib.exo.cdf_to_nc(args.cdfname)
