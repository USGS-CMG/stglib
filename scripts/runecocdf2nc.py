#!/usr/bin/env python

import stglib

args = stglib.cmd.ecocdf2nc_parser().parse_args()

ds = stglib.eco.cdf_to_nc(args.cdfname)
