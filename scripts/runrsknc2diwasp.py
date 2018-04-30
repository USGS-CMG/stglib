#!/usr/bin/env python

import stglib

args = stglib.cmd.rsknc2diwasp_parser().parse_args()

ds = stglib.rsk.nc2diwasp.nc_to_diwasp(args.ncname)
