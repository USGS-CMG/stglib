import argparse
import stglib
import xarray as xr


def cfcheckername_parser():
    description = "Run cfchecker on .nc files"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("ncname", help=".nc filename")

    return parser


args = cfcheckername_parser().parse_args()

print(args.ncname)

ds = xr.open_dataset(args.ncname)

stglib.core.utils.check_compliance(args.ncname, conventions=ds.attrs["Conventions"])
