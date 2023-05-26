import yaml

import stglib


def get_metadata(args):
    # initialize metadata from the globalatts file
    metadata = stglib.read_globalatts(args.gatts)

    # Add additional metadata from metadata config file
    with open(args.config) as f:
        config = yaml.safe_load(f)

    for k in config:
        metadata[k] = config[k]

    return metadata


def runaqdcdf2nc():
    args = stglib.cmd.aqdcdf2nc_parser().parse_args()

    if args.atmpres:
        ds = stglib.aqd.cdf2nc.cdf_to_nc(args.cdfname, atmpres=args.atmpres)
    else:
        ds = stglib.aqd.cdf2nc.cdf_to_nc(args.cdfname)


def runaqdhdr2cdf():
    args = stglib.cmd.aqdhdr2cdf_parser().parse_args()

    metadata = get_metadata(args)

    RAW = stglib.aqd.hdr2cdf.hdr_to_cdf(metadata)


def runaqdhrcdf2nc():
    args = stglib.cmd.aqdcdf2nc_parser().parse_args()

    if args.atmpres:
        ds = stglib.aqd.hrcdf2nc.cdf_to_nc(args.cdfname, atmpres=args.atmpres)
    else:
        ds = stglib.aqd.hrcdf2nc.cdf_to_nc(args.cdfname)


def runaqdhrhdr2cdf():
    args = stglib.cmd.aqdhdr2cdf_parser().parse_args()

    metadata = get_metadata(args)

    RAW = stglib.aqd.hrhdr2cdf.hdr_to_cdf(metadata)


def runecocdf2nc():
    args = stglib.cmd.ecocdf2nc_parser().parse_args()

    ds = stglib.eco.cdf_to_nc(args.cdfname)


def runecocsv2cdf():
    args = stglib.cmd.ecolog2cdf_parser().parse_args()

    metadata = get_metadata(args)

    RAW = stglib.eco.csv_to_cdf(metadata)


def runeofecdf2nc():
    args = stglib.cmd.eofecdf2nc_parser().parse_args()

    ds = stglib.eofe.cdf_to_nc(args.cdfname)


def runeofelog2cdf():
    args = stglib.cmd.eofelog2cdf_parser().parse_args()

    metadata = get_metadata(args)

    RAW = stglib.eofe.log_to_cdf(metadata)


def runexocdf2nc():
    args = stglib.cmd.exocdf2nc_parser().parse_args()

    if args.atmpres:
        ds = stglib.exo.cdf_to_nc(args.cdfname, atmpres=args.atmpres)
    else:
        ds = stglib.exo.cdf_to_nc(args.cdfname)


def runexocsv2cdf():
    args = stglib.cmd.exocsv2cdf_parser().parse_args()

    metadata = get_metadata(args)

    RAW = stglib.exo.csv_to_cdf(metadata)


def runhwlbcdf2nc():
    args = stglib.cmd.hwlbcdf2nc_parser().parse_args()

    ds = stglib.hobo.cdf_to_nc(args.cdfname)


def runhwlbcsv2cdf():
    args = stglib.cmd.hwlbcsv2cdf_parser().parse_args()

    metadata = get_metadata(args)

    RAW = stglib.hobo.csv_to_cdf(metadata)


def runiqcdf2nc():
    args = stglib.cmd.iqcdf2nc_parser().parse_args()

    ds = stglib.iq.cdf_to_nc(args.cdfname)


def runiqmat2cdf():
    args = stglib.cmd.iqmat2cdf_parser().parse_args()

    metadata = get_metadata(args)

    RAW = stglib.iq.mat_to_cdf(metadata)


def runrdicdf2nc():
    args = stglib.cmd.rdicdf2nc_parser().parse_args()

    if args.atmpres:
        ds = stglib.rdi.cdf2nc.cdf_to_nc(args.cdfname, atmpres=args.atmpres)
    else:
        ds = stglib.rdi.cdf2nc.cdf_to_nc(args.cdfname)


def runrdiraw2cdf():
    args = stglib.cmd.rdiraw2cdf_parser().parse_args()

    metadata = get_metadata(args)

    RAW = stglib.rdi.raw2cdf.raw_to_cdf(metadata)


def runrskcdf2nc():
    args = stglib.cmd.rskcdf2nc_parser().parse_args()

    if args.atmpres:
        ds = stglib.rsk.cdf2nc.cdf_to_nc(args.cdfname, atmpres=args.atmpres)
    else:
        ds = stglib.rsk.cdf2nc.cdf_to_nc(args.cdfname)


def runrskcsv2cdf():
    args = stglib.cmd.rskcsv2cdf_parser().parse_args()

    metadata = get_metadata(args)

    RAW = stglib.rsk.csv2cdf.csv_to_cdf(metadata)


def runrsknc2diwasp():
    args = stglib.cmd.rsknc2diwasp_parser().parse_args()

    ds = stglib.rsk.nc2diwasp.nc_to_diwasp(args.ncname)


def runrsknc2waves():
    args = stglib.cmd.rsknc2waves_parser().parse_args()

    ds = stglib.rsk.nc2waves.nc_to_waves(args.ncname)


def runrskrsk2cdf():
    args = stglib.cmd.rskrsk2cdf_parser().parse_args()

    metadata = get_metadata(args)

    RAW = stglib.rsk.rsk2cdf.rsk_to_cdf(metadata)


def runsigcdf2nc():
    args = stglib.cmd.sigcdf2nc_parser().parse_args()

    for f in args.cdfname:
        if args.atmpres:
            ds = stglib.sig.cdf2nc.cdf_to_nc(f, atmpres=args.atmpres)
        else:
            ds = stglib.sig.cdf2nc.cdf_to_nc(f)


def runsigdlfncdf2nc():
    args = stglib.cmd.sigdlfncdf2nc_parser().parse_args()

    for f in args.cdfname:
        if args.atmpres:
            ds = stglib.sig.dlfncdf2nc.cdf_to_nc(f, atmpres=args.atmpres)
        else:
            ds = stglib.sig.dlfncdf2nc.cdf_to_nc(f)


def runsigmat2cdf():
    args = stglib.cmd.sigmat2cdf_parser().parse_args()

    metadata = get_metadata(args)

    RAW = stglib.sig.mat2cdf.mat_to_cdf(metadata)


def runsigraw2cdf():
    args = stglib.cmd.sigraw2cdf_parser().parse_args()

    metadata = get_metadata(args)

    RAW = stglib.sig.raw2cdf.raw_to_cdf(metadata)


def runtrollcdf2nc():
    args = stglib.cmd.trollcdf2nc_parser().parse_args()

    ds = stglib.troll.cdf_to_nc(args.cdfname)


def runtrollcsv2cdf():
    args = stglib.cmd.trollcsv2cdf_parser().parse_args()

    metadata = get_metadata(args)

    RAW = stglib.troll.csv_to_cdf(metadata)


def runveccdf2nc():
    args = stglib.cmd.veccdf2nc_parser().parse_args()

    if args.atmpres:
        ds = stglib.vec.cdf2nc.cdf_to_nc(args.cdfname, atmpres=args.atmpres)
    else:
        ds = stglib.vec.cdf2nc.cdf_to_nc(args.cdfname)


def runvecdat2cdf():
    args = stglib.cmd.vechdr2cdf_parser().parse_args()

    metadata = get_metadata(args)

    RAW = stglib.vec.dat2cdf.dat_to_cdf(metadata)


def runwvscdf2nc():
    args = stglib.cmd.wvscdf2nc_parser().parse_args()

    if args.atmpres:
        VEL = stglib.aqd.wvscdf2nc.cdf_to_nc(args.cdfname, atmpres=args.atmpres)
    else:
        VEL = stglib.aqd.wvscdf2nc.cdf_to_nc(args.cdfname)


def runwvsnc2diwasp():
    args = stglib.cmd.wvsnc2diwasp_parser().parse_args()

    ds = stglib.aqd.wvsnc2diwasp.nc_to_diwasp(args.ncname)


def runwvsnc2waves():
    args = stglib.cmd.wvsnc2waves_parser().parse_args()

    ds = stglib.aqd.wvsnc2waves.nc_to_waves(args.ncname)


def runwvswad2cdf():
    args = stglib.cmd.wvswad2cdf_parser().parse_args()

    metadata = get_metadata(args)

    RAW = stglib.aqd.wvswad2cdf.wad_to_cdf(metadata)


def runwxtcdf2nc():
    import stglib.wxt

    args = stglib.cmd.wxtcdf2nc_parser().parse_args()

    ds = stglib.wxt.cdf_to_nc(args.cdfname)


def runwxtcsv2cdf():
    import stglib.wxt

    args = stglib.cmd.wxtcsv2cdf_parser().parse_args()

    metadata = get_metadata(args)

    RAW = stglib.wxt.csv_to_cdf(metadata)
