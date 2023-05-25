import yaml

import stglib


def runsigmat2cdf():
    args = stglib.cmd.sigmat2cdf_parser().parse_args()

    # initialize metadata from the globalatts file
    metadata = stglib.read_globalatts(args.gatts)

    # Add additional metadata from metadata config file
    with open(args.config) as f:
        config = yaml.safe_load(f)

    for k in config:
        metadata[k] = config[k]

    RAW = stglib.sig.mat2cdf.mat_to_cdf(metadata)


def runsigcdf2nc():
    args = stglib.cmd.sigcdf2nc_parser().parse_args()

    for f in args.cdfname:
        if args.atmpres:
            ds = stglib.sig.cdf2nc.cdf_to_nc(f, atmpres=args.atmpres)
        else:
            ds = stglib.sig.cdf2nc.cdf_to_nc(f)


def runvecdat2cdf():
    args = stglib.cmd.vechdr2cdf_parser().parse_args()

    # initialize metadata from the globalatts file
    metadata = stglib.read_globalatts(args.gatts)

    # Add additional metadata from metadata config file
    with open(args.config) as f:
        config = yaml.safe_load(f)

    for k in config:
        metadata[k] = config[k]

    RAW = stglib.vec.dat2cdf.dat_to_cdf(metadata)


def runveccdf2nc():
    args = stglib.cmd.veccdf2nc_parser().parse_args()

    if args.atmpres:
        ds = stglib.vec.cdf2nc.cdf_to_nc(args.cdfname, atmpres=args.atmpres)
    else:
        ds = stglib.vec.cdf2nc.cdf_to_nc(args.cdfname)
