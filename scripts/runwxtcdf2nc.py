#!/usr/bin/env python

import stglib
import stglib.wxt

args = stglib.cmd.wxtcdf2nc_parser().parse_args()

ds = stglib.wxt.cdf_to_nc(args.cdfname)
