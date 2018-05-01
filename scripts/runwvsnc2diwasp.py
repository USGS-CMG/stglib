#!/usr/bin/env python

import stglib

args = stglib.cmd.wvsnc2diwasp_parser().parse_args()

ds = stglib.aqd.wvsnc2diwasp.nc_to_diwasp(args.ncname)
