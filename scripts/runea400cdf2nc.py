#!/usr/bin/env python

import stglib

args = stglib.cmd.ea400cdf2nc_parser().parse_args()

ds = stglib.ea400.cdf_to_nc(args.cdfname)
