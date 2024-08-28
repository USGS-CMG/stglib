import warnings

import yaml

import stglib


def get_metadata(args):
    # initialize metadata from the globalatts file
    metadata = stglib.read_globalatts(args.gatts)

    # Add additional metadata from metadata config file
    with open(args.config) as f:
        config = yaml.safe_load(f)

    try:
        for k in config:
            if k in metadata:
                warnings.warn(
                    f"attrs collision. Replacing '{k}={metadata[k]}' from global attributes file with '{k}={config[k]}' from YAML config file."
                )
            metadata[k] = config[k]
    except TypeError:
        raise TypeError(f"Could not load metadata from {args.config}")

    return metadata


def run_cdf_to_nc(f, args):
    # if not a list, make it one so we can process multiple input files if necessary
    if not isinstance(args.cdfname, list):
        args.cdfname = [args.cdfname]

    for cdfname in args.cdfname:
        if hasattr(args, "atmpres") and args.atmpres:
            ds = f(cdfname, atmpres=args.atmpres)
        else:
            ds = f(cdfname)

    return ds


def runaqdcdf2nc(args=None):
    if not args:
        args = stglib.cmd.aqdcdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.aqd.cdf2nc.cdf_to_nc, args)


def runaqdhdr2cdf(args=None):
    if not args:
        args = stglib.cmd.aqdhdr2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.aqd.hdr2cdf.hdr_to_cdf(metadata)


def runaqdhrcdf2nc(args=None):
    if not args:
        args = stglib.cmd.aqdcdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.aqd.hrcdf2nc.cdf_to_nc, args)


def runaqdhrhdr2cdf(args=None):
    if not args:
        args = stglib.cmd.aqdhdr2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.aqd.hrhdr2cdf.hdr_to_cdf(metadata)


def runecocdf2nc(args=None):
    if not args:
        args = stglib.cmd.ecocdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.eco.cdf_to_nc, args)


def runecocsv2cdf(args=None):
    if not args:
        args = stglib.cmd.ecolog2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.eco.csv_to_cdf(metadata)


def runeofecdf2nc(args=None):
    if not args:
        args = stglib.cmd.eofecdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.eofe.cdf_to_nc, args)


def runeofelog2cdf(args=None):
    if not args:
        args = stglib.cmd.eofelog2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.eofe.log_to_cdf(metadata)


def runexocdf2nc(args=None):
    if not args:
        args = stglib.cmd.exocdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.exo.cdf_to_nc, args)


def runexocsv2cdf(args=None):
    if not args:
        args = stglib.cmd.exocsv2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.exo.csv_to_cdf(metadata)


def runhobocdf2nc(args=None):
    if not args:
        args = stglib.cmd.hobocdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.hobo.cdf_to_nc, args)


def runhobocsv2cdf(args=None):
    if not args:
        args = stglib.cmd.hobocsv2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.hobo.csv_to_cdf(metadata)


def runiqcdf2nc(args=None):
    if not args:
        args = stglib.cmd.iqcdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.iq.cdf_to_nc, args)


def runiqmat2cdf(args=None):
    if not args:
        args = stglib.cmd.iqmat2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.iq.mat_to_cdf(metadata)


def runlisstcdf2nc(args=None):
    if not args:
        args = stglib.cmd.lisstcdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.lisst.cdf_to_nc, args)


def runlisstcsv2cdf(args=None):
    if not args:
        args = stglib.cmd.lisstcsv2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.lisst.csv_to_cdf(metadata)


def runrdicdf2nc(args=None):
    if not args:
        args = stglib.cmd.rdicdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.rdi.cdf2nc.cdf_to_nc, args)


def runrdiraw2cdf(args=None):
    if not args:
        args = stglib.cmd.rdiraw2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.rdi.raw2cdf.raw_to_cdf(metadata)


def runrskcdf2nc(args=None):
    if not args:
        args = stglib.cmd.rskcdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.rsk.cdf2nc.cdf_to_nc, args)


def runrskcsv2cdf(args=None):
    if not args:
        args = stglib.cmd.rskcsv2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.rsk.csv2cdf.csv_to_cdf(metadata)


def runrsknc2waves(args=None):
    if not args:
        args = stglib.cmd.rsknc2waves_parser().parse_args()

    stglib.rsk.nc2waves.nc_to_waves(args.ncname)


def runrskrsk2cdf(args=None):
    if not args:
        args = stglib.cmd.rskrsk2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.rsk.rsk2cdf.rsk_to_cdf(metadata)


def runsigcdf2nc(args=None):
    if not args:
        args = stglib.cmd.sigcdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.sig.cdf2nc.cdf_to_nc, args)


def runsigmat2cdf(args=None):
    if not args:
        args = stglib.cmd.sigmat2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.sig.mat2cdf.mat_to_cdf(metadata)


def runtrollcdf2nc(args=None):
    if not args:
        args = stglib.cmd.trollcdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.troll.cdf_to_nc, args)


def runtrollcsv2cdf(args=None):
    if not args:
        args = stglib.cmd.trollcsv2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.troll.csv_to_cdf(metadata)


def runveccdf2nc(args=None):
    if not args:
        args = stglib.cmd.veccdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.vec.cdf2nc.cdf_to_nc, args)


def runvecdat2cdf(args=None):
    if not args:
        args = stglib.cmd.vechdr2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.vec.dat2cdf.dat_to_cdf(metadata)


def runvecnc2waves(args=None):
    if not args:
        args = stglib.cmd.vecnc2waves_parser().parse_args()

    stglib.vec.nc2waves.nc_to_waves(args.ncname)


def runwvscdf2nc(args=None):
    if not args:
        args = stglib.cmd.wvscdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.aqd.wvscdf2nc.cdf_to_nc, args)


def runwvsnc2diwasp(args=None):
    if not args:
        args = stglib.cmd.wvsnc2diwasp_parser().parse_args()

    stglib.aqd.wvsnc2diwasp.nc_to_diwasp(args.ncname)


def runwvsnc2waves(args=None):
    if not args:
        args = stglib.cmd.wvsnc2waves_parser().parse_args()

    stglib.aqd.wvsnc2waves.nc_to_waves(args.ncname)


def runwvswad2cdf(args=None):
    if not args:
        args = stglib.cmd.wvswad2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.aqd.wvswad2cdf.wad_to_cdf(metadata)


def runwxtcdf2nc(args=None):
    if not args:
        args = stglib.cmd.wxtcdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.wxt.cdf_to_nc, args)


def runwxtcsv2cdf(args=None):
    if not args:
        args = stglib.cmd.wxtcsv2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.wxt.csv_to_cdf(metadata)


def runtcmcdf2nc(args=None):
    if not args:
        args = stglib.cmd.tcmcdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.tcm.cdf_to_nc, args)


def runtcmcsv2cdf(args=None):
    if not args:
        args = stglib.cmd.tcmcsv2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.tcm.csv_to_cdf(metadata)


def runmccdf2nc(args=None):
    if not args:
        args = stglib.cmd.mccdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.mc.cdf_to_nc, args)


def runmcasc2cdf(args=None):
    if not args:
        args = stglib.cmd.mcasc2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.mc.asc_to_cdf(metadata)


def runsgcdf2nc(args=None):
    if not args:
        args = stglib.cmd.sgcdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.sg.cdf_to_nc, args)


def runsgtid2cdf(args=None):
    if not args:
        args = stglib.cmd.sgtid2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.sg.tid_to_cdf(metadata)


def runsgnc2waves(args=None):
    if not args:
        args = stglib.cmd.sgnc2waves_parser().parse_args()

    stglib.sg.nc_to_waves(args.ncname)


def runots():
    args = stglib.cmd.runots_parser().parse_args()

    if args.instrument == "aqd":
        if args.step == "hdr2cdf":
            runaqdhdr2cdf(args)
        elif args.step == "cdf2nc":
            runaqdcdf2nc(args)
    if args.instrument == "aqdhr":
        if args.step == "hdr2cdf":
            runaqdhrhdr2cdf(args)
        elif args.step == "cdf2nc":
            runaqdhrcdf2nc(args)
    elif args.instrument in ["aqdwvs", "wvs"]:
        if args.step == "wad2cdf":
            runwvsdat2cdf(args)
        elif args.step == "cdf2nc":
            runwvscdf2nc(args)
        elif args.step == "nc2waves":
            runwvsnc2waves(args)
    elif args.instrument in ["rbr", "rsk"]:
        if args.step == "csv2cdf":
            runrskcsv2cdf(args)
        elif args.step == "cdf2nc":
            runrskcdf2nc(args)
        elif args.step == "nc2waves":
            runrsknc2waves(args)
    elif args.instrument == "sig":
        if args.step == "mat2cdf":
            runsigmat2cdf(args)
        elif args.step == "cdf2nc":
            runsigcdf2nc(args)
    elif args.instrument == "vec":
        if args.step == "dat2cdf":
            runvecdat2cdf(args)
        elif args.step == "cdf2nc":
            runveccdf2nc(args)
        elif args.step == "nc2waves":
            runvecnc2waves(args)
    elif args.instrument == "eco":
        if args.step == "csv2cdf":
            runecocsv2cdf(args)
        elif args.step == "cdf2nc":
            runecocdf2nc(args)
    elif args.instrument == "eofe":
        if args.step == "log2cdf":
            runeofelog2cdf(args)
        elif args.step == "cdf2nc":
            runeofecdf2nc(args)
    elif args.instrument == "exo":
        if args.step == "csv2cdf":
            runexomat2cdf(args)
        elif args.step == "cdf2nc":
            runexocdf2nc(args)
    elif args.instrument == "hobo":
        if args.step == "csv2cdf":
            runhobocsv2cdf(args)
        elif args.step == "cdf2nc":
            runhobocdf2nc(args)
    elif args.instrument == "iq":
        if args.step == "mat2cdf":
            runiqmat2cdf(args)
        elif args.step == "cdf2nc":
            runiqcdf2nc(args)
    elif args.instrument == "lisst":
        if args.step == "csv2cdf":
            runlisstcsv2cdf(args)
        elif args.step == "cdf2nc":
            runlisstcdf2nc(args)
    elif args.instrument == "mc":
        if args.step == "asc2cdf":
            runmcasc2cdf(args)
        elif args.step == "cdf2nc":
            runmccdf2nc(args)
    elif args.instrument == "sg":
        if args.step == "tid2cdf":
            runsgtid2cdf(args)
        elif args.step == "cdf2nc":
            runsgcdf2nc(args)
        elif args.step == "nc2waves":
            runsgnc2waves(args)
    elif args.instrument == "tcm":
        if args.step == "csv2cdf":
            runtcmcsv2cdf(args)
        elif args.step == "cdf2nc":
            runtcmcdf2nc(args)
    elif args.instrument == "troll":
        if args.step == "csv2cdf":
            runtrollcsv2cdf(args)
        elif args.step == "cdf2nc":
            runtrollcdf2nc(args)
    elif args.instrument == "wxt":
        if args.step == "csv2cdf":
            runwxtcsv2cdf(args)
        elif args.step == "cdf2nc":
            runwxtcdf2nc(args)
