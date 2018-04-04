#!/usr/bin/env python

import sys
sys.path.insert(0, '/Users/dnowacki/Documents/stglib')
import stglib
import argparse

parser = argparse.ArgumentParser(description='Convert raw Aquadopp .cdf wave files to processed .nc files')
parser.add_argument('cdfname', help='raw .CDF filename')
parser.add_argument('--atmpres', help='path to cdf file containing atmopsheric pressure data')

args = parser.parse_args()

if args.atmpres:
    VEL = stglib.aqd.wvscdf2nc.cdf_to_nc(args.cdfname, atmpres=args.atmpres)
else:
    VEL = stglib.aqd.wvscdf2nc.cdf_to_nc(args.cdfname)
