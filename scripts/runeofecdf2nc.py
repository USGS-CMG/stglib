#!/usr/bin/env python

import stglib

args = stglib.cmd.eofecdf2nc_parser().parse_args()

ds = stglib.eofe.cdf_to_nc(args.cdfname)
