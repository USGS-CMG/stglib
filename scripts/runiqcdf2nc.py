#!/usr/bin/env python

import stglib

args = stglib.cmd.iqcdf2nc_parser().parse_args()

ds = stglib.iq.cdf_to_nc(args.cdfname)
