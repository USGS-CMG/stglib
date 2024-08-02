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


def runaqdcdf2nc():
    args = stglib.cmd.aqdcdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.aqd.cdf2nc.cdf_to_nc, args)


def runaqdhdr2cdf():
    args = stglib.cmd.aqdhdr2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.aqd.hdr2cdf.hdr_to_cdf(metadata)


def runaqdhrcdf2nc():
    args = stglib.cmd.aqdcdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.aqd.hrcdf2nc.cdf_to_nc, args)


def runaqdhrhdr2cdf():
    args = stglib.cmd.aqdhdr2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.aqd.hrhdr2cdf.hdr_to_cdf(metadata)


def runecocdf2nc():
    args = stglib.cmd.ecocdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.eco.cdf_to_nc, args)


def runecocsv2cdf():
    args = stglib.cmd.ecolog2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.eco.csv_to_cdf(metadata)


def runeofecdf2nc():
    args = stglib.cmd.eofecdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.eofe.cdf_to_nc, args)


def runeofelog2cdf():
    args = stglib.cmd.eofelog2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.eofe.log_to_cdf(metadata)


def runexocdf2nc():
    args = stglib.cmd.exocdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.exo.cdf_to_nc, args)


def runexocsv2cdf():
    args = stglib.cmd.exocsv2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.exo.csv_to_cdf(metadata)


def runhobocdf2nc():
    args = stglib.cmd.hobocdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.hobo.cdf_to_nc, args)


def runhobocsv2cdf():
    args = stglib.cmd.hobocsv2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.hobo.csv_to_cdf(metadata)


def runiqcdf2nc():
    args = stglib.cmd.iqcdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.iq.cdf_to_nc, args)


def runiqmat2cdf():
    args = stglib.cmd.iqmat2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.iq.mat_to_cdf(metadata)


def runlisstcdf2nc():
    args = stglib.cmd.lisstcdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.lisst.cdf_to_nc, args)


def runlisstcsv2cdf():
    args = stglib.cmd.lisstcsv2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.lisst.csv_to_cdf(metadata)


def runrdicdf2nc():
    args = stglib.cmd.rdicdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.rdi.cdf2nc.cdf_to_nc, args)


def runrdiraw2cdf():
    args = stglib.cmd.rdiraw2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.rdi.raw2cdf.raw_to_cdf(metadata)


def runrskcdf2nc():
    args = stglib.cmd.rskcdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.rsk.cdf2nc.cdf_to_nc, args)


def runrskcsv2cdf():
    args = stglib.cmd.rskcsv2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.rsk.csv2cdf.csv_to_cdf(metadata)


def runrsknc2waves():
    args = stglib.cmd.rsknc2waves_parser().parse_args()

    stglib.rsk.nc2waves.nc_to_waves(args.ncname)


def runrskrsk2cdf():
    args = stglib.cmd.rskrsk2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.rsk.rsk2cdf.rsk_to_cdf(metadata)


def runsigcdf2nc():
    args = stglib.cmd.sigcdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.sig.cdf2nc.cdf_to_nc, args)


def runsigmat2cdf():
    args = stglib.cmd.sigmat2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.sig.mat2cdf.mat_to_cdf(metadata)


def runtrollcdf2nc():
    args = stglib.cmd.trollcdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.troll.cdf_to_nc, args)


def runtrollcsv2cdf():
    args = stglib.cmd.trollcsv2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.troll.csv_to_cdf(metadata)


def runveccdf2nc():
    args = stglib.cmd.veccdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.vec.cdf2nc.cdf_to_nc, args)


def runvecdat2cdf():
    args = stglib.cmd.vechdr2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.vec.dat2cdf.dat_to_cdf(metadata)


def runvecnc2waves():
    args = stglib.cmd.vecnc2waves_parser().parse_args()

    stglib.vec.nc2waves.nc_to_waves(args.ncname)


def runwvscdf2nc():
    args = stglib.cmd.wvscdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.aqd.wvscdf2nc.cdf_to_nc, args)


def runwvsnc2diwasp():
    args = stglib.cmd.wvsnc2diwasp_parser().parse_args()

    stglib.aqd.wvsnc2diwasp.nc_to_diwasp(args.ncname)


def runwvsnc2waves():
    args = stglib.cmd.wvsnc2waves_parser().parse_args()

    stglib.aqd.wvsnc2waves.nc_to_waves(args.ncname)


def runwvswad2cdf():
    args = stglib.cmd.wvswad2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.aqd.wvswad2cdf.wad_to_cdf(metadata)


def runwxtcdf2nc():
    args = stglib.cmd.wxtcdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.wxt.cdf_to_nc, args)


def runwxtcsv2cdf():
    args = stglib.cmd.wxtcsv2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.wxt.csv_to_cdf(metadata)


def runtcmcdf2nc():
    args = stglib.cmd.tcmcdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.tcm.cdf_to_nc, args)


def runtcmcsv2cdf():
    args = stglib.cmd.tcmcsv2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.tcm.csv_to_cdf(metadata)


def runmccdf2nc():
    args = stglib.cmd.mccdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.mc.cdf_to_nc, args)


def runmcasc2cdf():
    args = stglib.cmd.mcasc2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.mc.asc_to_cdf(metadata)


def runsgcdf2nc():
    args = stglib.cmd.sgcdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.sg.cdf_to_nc, args)


def runsgtid2cdf():
    args = stglib.cmd.sgtid2cdf_parser().parse_args()

    metadata = get_metadata(args)

    stglib.sg.tid_to_cdf(metadata)


def runsgnc2waves():
    args = stglib.cmd.sgnc2waves_parser().parse_args()

    stglib.sg.nc_to_waves(args.ncname)
