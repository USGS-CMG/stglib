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
    except TypeError as e:
        e.add_note(f"Could not load metadata from {args.config}")
        raise

    return metadata


def run_cdf_to_nc(f, args):
    # if not a list, make it one so we can process multiple input files if necessary
    if not isinstance(args.cdfname, list):
        args.cdfname = [args.cdfname]

    kwargs = {
        k: getattr(args, k)
        for k in ("atmpres", "salwtemp", "height")
        if getattr(args, k, None)
    }
    for cdfname in args.cdfname:
        ds = f(cdfname, **kwargs)
    return ds


def runabssmat2cdf(args=None):
    if not args:
        args = stglib.cmd.inst2cdf_parser(
            "Convert ABS .mat files to raw .cdf format. Run this script from the directory containing ABS glob_att and config files."
        ).parse_args()

    metadata = get_metadata(args)

    stglib.abss.mat2cdf(metadata)


def runabsscdf2nc(args=None):
    if not args:
        args = stglib.cmd.cdf2nc_parser(
            "Convert raw .cdf format to processed .nc files, optionally compensating for atmospheric pressure"
        ).parse_args()

    run_cdf_to_nc(stglib.abss.cdf2nc, args)


def runaqdcdf2nc(args=None):
    if not args:
        args = stglib.cmd.cdf2nc_parser(
            "Convert raw Aquadopp .cdf format to processed .nc files, optionally compensating for atmospheric pressure"
        ).parse_args()

    run_cdf_to_nc(stglib.aqd.cdf2nc.cdf_to_nc, args)


def runaqdhdr2cdf(args=None):
    if not args:
        args = stglib.cmd.inst2cdf_parser(
            "Convert Aquadopp text files to raw .cdf format. Run this script from the directory containing Aquadopp files."
        ).parse_args()

    metadata = get_metadata(args)

    stglib.aqd.hdr2cdf.hdr_to_cdf(metadata)


def runaqdhrcdf2nc(args=None):
    if not args:
        args = stglib.cmd.cdf2nc_parser(
            "Convert raw Aquadopp HR .cdf format to processed .nc files, optionally compensating for atmospheric pressure"
        ).parse_args()

    run_cdf_to_nc(stglib.aqd.hrcdf2nc.cdf_to_nc, args)


def runaqdhrhdr2cdf(args=None):
    if not args:
        args = stglib.cmd.inst2cdf_parser(
            "Convert Aquadopp HR text files to raw .cdf format. Run this script from the directory containing Aquadopp files."
        ).parse_args()

    metadata = get_metadata(args)

    stglib.aqd.hrhdr2cdf.hdr_to_cdf(metadata)


def runecocdf2nc(args=None):
    if not args:
        args = stglib.cmd.cdf2nc_parser(
            "Convert raw ECO .cdf format to processed .nc files", atmpres=False
        ).parse_args()

    run_cdf_to_nc(stglib.eco.cdf_to_nc, args)


def runecocsv2cdf(args=None):
    if not args:
        args = stglib.cmd.inst2cdf_parser(
            "Convert WET Labs ECO file to raw .cdf format. Run this script from the directory containing ECO file."
        ).parse_args()

    metadata = get_metadata(args)

    stglib.eco.csv_to_cdf(metadata)


def runeofecdf2nc(args=None):
    if not args:
        args = stglib.cmd.cdf2nc_parser(
            "Convert raw echologger .cdf format to processed .nc files", atmpres=False
        ).parse_args()

    run_cdf_to_nc(stglib.eofe.cdf_to_nc, args)


def runeofelog2cdf(args=None):
    if not args:
        args = stglib.cmd.inst2cdf_parser(
            "Convert EofE echologger .log file to raw .cdf format. Run this script from the directory containing ea400 echologger .log file."
        ).parse_args()

    metadata = get_metadata(args)

    stglib.eofe.log_to_cdf(metadata)


def runexocdf2nc(args=None):
    if not args:
        args = stglib.cmd.cdf2nc_parser(
            "Convert raw EXO .cdf format to processed .nc files"
        ).parse_args()

    run_cdf_to_nc(stglib.exo.cdf_to_nc, args)


def runexocsv2cdf(args=None):
    if not args:
        args = stglib.cmd.inst2cdf_parser(
            "Convert EXO .csv file to raw .cdf format. Run this script from the directory containing EXO file."
        ).parse_args()

    metadata = get_metadata(args)

    stglib.exo.csv_to_cdf(metadata)


def runglxcdf2nc(args=None):
    if not args:
        args = stglib.cmd.cdf2nc_parser(
            "Convert raw Geolux wave radar .cdf format to processed .nc files",
            atmpres=False,
        ).parse_args()

    run_cdf_to_nc(stglib.glx.cdf_to_nc, args)


def runglxdat2cdf(args=None):
    if not args:
        args = stglib.cmd.inst2cdf_parser(
            "Convert exported Geolux wave radar files to raw .cdf format. Run this script from the directory containing the files."
        ).parse_args()

    metadata = get_metadata(args)

    stglib.glx.dat_to_cdf(metadata)


def runglxnc2waves(args=None):
    if not args:
        args = stglib.cmd.nc2waves_parser("Generate waves statistics file").parse_args()

    stglib.glx.nc_to_waves(args.ncname)


def runhobocdf2nc(args=None):
    if not args:
        args = stglib.cmd.cdf2nc_parser(
            "Convert raw HOBO .cdf format to processed .nc files", atmpres=False
        ).parse_args()

    run_cdf_to_nc(stglib.hobo.cdf_to_nc, args)


def runhobocsv2cdf(args=None):
    if not args:
        args = stglib.cmd.inst2cdf_parser(
            "Convert HOBO pressure sensor .csv file to raw .cdf format. Run this script from the directory containing HOBO file."
        ).parse_args()

    metadata = get_metadata(args)

    stglib.hobo.csv_to_cdf(metadata)


def runiqcdf2nc(args=None):
    if not args:
        args = stglib.cmd.cdf2nc_parser(
            "Convert raw SonTek IQ .cdf format to processed .nc files", atmpres=False
        ).parse_args()

    run_cdf_to_nc(stglib.iq.cdf_to_nc, args)


def runiqmat2cdf(args=None):
    if not args:
        args = stglib.cmd.inst2cdf_parser(
            "Convert SonTek IQ .mat file to raw .cdf format. Run this script from the directory containing IQ file."
        ).parse_args()

    metadata = get_metadata(args)

    stglib.iq.mat_to_cdf(metadata)


def runlisstcdf2nc(args=None):
    if not args:
        args = stglib.cmd.cdf2nc_parser(
            "Convert raw LISST .cdf format to processed .nc files"
        ).parse_args()

    run_cdf_to_nc(stglib.lisst.cdf_to_nc, args)


def runlisstcsv2cdf(args=None):
    if not args:
        args = stglib.cmd.inst2cdf_parser(
            "Convert LISST .csv file to raw .cdf format. Run this script from the directory containing LISST file."
        ).parse_args()

    metadata = get_metadata(args)

    stglib.lisst.csv_to_cdf(metadata)


def runrdicdf2nc(args=None):
    if not args:
        args = stglib.cmd.cdf2nc_parser(
            "Convert raw RDI .cdf format to processed .nc files, optionally compensating for atmospheric pressure"
        ).parse_args()

    run_cdf_to_nc(stglib.rdi.cdf2nc.cdf_to_nc, args)


def runrdiraw2cdf(args=None):
    if not args:
        args = stglib.cmd.inst2cdf_parser(
            "Convert RDI raw binary files to raw .cdf format. Run this script from the directory containing RDI files."
        ).parse_args()

    metadata = get_metadata(args)

    stglib.rdi.raw2cdf.raw_to_cdf(metadata)


def runrskcdf2nc(args=None):
    if not args:
        args = stglib.cmd.cdf2nc_parser(
            "Convert raw RBR d|wave .cdf format to processed .nc files, optionally compensating for atmospheric pressure"
        ).parse_args()

    run_cdf_to_nc(stglib.rsk.cdf2nc.cdf_to_nc, args)


def runrskcsv2cdf(args=None):
    if not args:
        args = stglib.cmd.inst2cdf_parser(
            "Convert exported RBR csv files to raw .cdf format. Run this script from the directory containing the files."
        ).parse_args()

    metadata = get_metadata(args)

    stglib.rsk.csv2cdf.csv_to_cdf(metadata)


def runrsknc2waves(args=None):
    if not args:
        args = stglib.cmd.nc2waves_parser("Generate waves statistics file").parse_args()

    if hasattr(args, "salwtemp") and args.salwtemp:
        stglib.rsk.nc2waves.nc_to_waves(args.ncname, salwtemp=args.salwtemp)
    else:
        stglib.rsk.nc2waves.nc_to_waves(args.ncname)


def runrsknc2diwasp(args=None):
    if not args:
        args = stglib.cmd.nc2diwasp_parser(
            "Generate DIWASP waves statistics file"
        ).parse_args()

    if hasattr(args, "salwtemp") and args.salwtemp:
        stglib.rsk.nc2waves.nc_to_diwasp(args.ncname, salwtemp=args.salwtemp)
    else:
        stglib.rsk.nc2waves.nc_to_diwasp(args.ncname)


def runrskrsk2cdf(args=None):
    if not args:
        args = stglib.cmd.inst2cdf_parser(
            "Convert raw RBR files (.rsk) to raw .cdf format. Run this script from the directory containing the files."
        ).parse_args()

    metadata = get_metadata(args)

    stglib.rsk.rsk2cdf.rsk_to_cdf(metadata)


def runsignc2diwasp(args=None):
    if not args:
        args = stglib.cmd.nc2diwasp_parser(
            "Generate DIWASP waves statistics file"
        ).parse_args()

    if hasattr(args, "salwtemp") and args.salwtemp:
        stglib.sig.nc2waves.nc_to_diwasp(args.ncname, salwtemp=args.salwtemp)
    else:
        stglib.sig.nc2waves.nc_to_diwasp(args.ncname)


def runsignc2waves(args=None):
    if not args:
        args = stglib.cmd.nc2waves_parser("Generate waves statistics file").parse_args()

    if hasattr(args, "salwtemp") and args.salwtemp:
        stglib.sig.nc2waves.nc_to_waves(args.ncname, salwtemp=args.salwtemp)
    else:
        stglib.sig.nc2waves.nc_to_waves(args.ncname)


def runsigcdf2nc(args=None):
    if not args:
        args = stglib.cmd.sigcdf2nc_parser().parse_args()

    run_cdf_to_nc(stglib.sig.cdf2nc.cdf_to_nc, args)


def runsigmat2cdf(args=None):
    if not args:
        args = stglib.cmd.inst2cdf_parser(
            "Convert Signature files exported in Matlab format to raw .cdf format. Run this script from the directory containing Signature files."
        ).parse_args()

    metadata = get_metadata(args)

    stglib.sig.mat2cdf.mat_to_cdf(metadata)


def runtrollcdf2nc(args=None):
    if not args:
        args = stglib.cmd.cdf2nc_parser(
            "Convert raw Aqua TROLL .cdf format to processed .nc files", atmpres=False
        ).parse_args()

    run_cdf_to_nc(stglib.troll.cdf_to_nc, args)


def runtrollcsv2cdf(args=None):
    if not args:
        args = stglib.cmd.inst2cdf_parser(
            "Convert Aqua TROLL .csv file to raw .cdf format. Run this script from the directory containing .csv file."
        ).parse_args()

    metadata = get_metadata(args)

    stglib.troll.csv_to_cdf(metadata)


def runveccdf2nc(args=None):
    if not args:
        args = stglib.cmd.cdf2nc_parser(
            "Convert raw Vector .cdf format to processed .nc files, optionally compensating for atmospheric pressure"
        ).parse_args()

    run_cdf_to_nc(stglib.vec.cdf2nc.cdf_to_nc, args)


def runvecdat2cdf(args=None):
    if not args:
        args = stglib.cmd.inst2cdf_parser(
            "Convert Vector text files to raw .cdf format. Run this script from the directory containing Vector files."
        ).parse_args()

    metadata = get_metadata(args)

    stglib.vec.dat2cdf.dat_to_cdf(metadata)


def runvecnc2waves(args=None):
    if not args:
        args = stglib.cmd.nc2waves_parser(
            "Generate Vector waves statistics file"
        ).parse_args()

    stglib.vec.nc2waves.nc_to_waves(args.ncname)


def runwvscdf2nc(args=None):
    if not args:
        args = stglib.cmd.cdf2nc_parser(
            "Convert raw Aquadopp .cdf wave files to processed .nc files"
        ).parse_args()

    run_cdf_to_nc(stglib.aqd.wvscdf2nc.cdf_to_nc, args)


def runwvsnc2diwasp(args=None):
    if not args:
        args = stglib.cmd.nc2diwasp_parser(
            "Convert processed Aquadopp waves .nc files using DIWASP"
        ).parse_args()

    stglib.aqd.wvsnc2diwasp.nc_to_diwasp(args.ncname)


def runwvsnc2waves(args=None):
    if not args:
        args = stglib.cmd.nc2waves_parser("Generate waves statistics file").parse_args()

    stglib.aqd.wvsnc2waves.nc_to_waves(args.ncname)


def runwvswad2cdf(args=None):
    if not args:
        args = stglib.cmd.inst2cdf_parser(
            "Convert Aquadopp .wad wave files to raw .cdf format. Run this script from the directory containing Aquadopp files."
        ).parse_args()

    metadata = get_metadata(args)

    stglib.aqd.wvswad2cdf.wad_to_cdf(metadata)


def runmetcdf2nc(args=None):
    if not args:
        args = stglib.cmd.cdf2nc_parser(
            "Convert raw met station .cdf format to processed .nc files", atmpres=False
        ).parse_args()

    run_cdf_to_nc(stglib.met.cdf_to_nc, args)


def runmetcsv2cdf(args=None):
    if not args:
        args = stglib.cmd.inst2cdf_parser(
            "Convert met .csv file to raw .cdf format. Run this script from the directory containing .csv file."
        ).parse_args()

    metadata = get_metadata(args)

    stglib.met.csv_to_cdf(metadata)


def runtcmcdf2nc(args=None):
    if not args:
        args = stglib.cmd.cdf2nc_parser(
            "Convert raw TCM .cdf format to processed .nc files", atmpres=False
        ).parse_args()

    run_cdf_to_nc(stglib.tcm.cdf_to_nc, args)


def runtcmcsv2cdf(args=None):
    if not args:
        args = stglib.cmd.inst2cdf_parser(
            "Convert Lowell Tilt Current Meter .txt file to raw .cdf format. Run this script from the directory containing TCM file."
        ).parse_args()

    metadata = get_metadata(args)

    stglib.tcm.csv_to_cdf(metadata)


def runmccdf2nc(args=None):
    if not args:
        args = stglib.cmd.cdf2nc_parser(
            "Convert raw SBE 37 MicroCAT .cdf format to processed .nc files",
            atmpres=False,
        ).parse_args()

    run_cdf_to_nc(stglib.mc.cdf_to_nc, args)


def runmcasc2cdf(args=None):
    if not args:
        args = stglib.cmd.inst2cdf_parser(
            "Convert SBE 37 MicroCAT .asc file to raw .cdf format. Run this script from the directory containing MicroCAT .asc file."
        ).parse_args()

    metadata = get_metadata(args)

    stglib.mc.asc_to_cdf(metadata)


_REGISTRY = {
    ("abss", "mat2cdf"): runabssmat2cdf,
    ("abss", "cdf2nc"): runabsscdf2nc,
    ("aqd", "hdr2cdf"): runaqdhdr2cdf,
    ("aqd", "cdf2nc"): runaqdcdf2nc,
    ("aqdhr", "hdr2cdf"): runaqdhrhdr2cdf,
    ("aqdhr", "cdf2nc"): runaqdhrcdf2nc,
    ("aqdwvs", "wad2cdf"): runwvswad2cdf,
    ("aqdwvs", "cdf2nc"): runwvscdf2nc,
    ("aqdwvs", "nc2waves"): runwvsnc2waves,
    ("wvs", "wad2cdf"): runwvswad2cdf,
    ("wvs", "cdf2nc"): runwvscdf2nc,
    ("wvs", "nc2waves"): runwvsnc2waves,
    ("rbr", "csv2cdf"): runrskcsv2cdf,
    ("rbr", "cdf2nc"): runrskcdf2nc,
    ("rbr", "nc2waves"): runrsknc2waves,
    ("rbr", "nc2diwasp"): runrsknc2diwasp,
    ("rsk", "csv2cdf"): runrskcsv2cdf,
    ("rsk", "cdf2nc"): runrskcdf2nc,
    ("rsk", "nc2waves"): runrsknc2waves,
    ("rsk", "nc2diwasp"): runrsknc2diwasp,
    ("sig", "mat2cdf"): runsigmat2cdf,
    ("sig", "cdf2nc"): runsigcdf2nc,
    ("sig", "nc2waves"): runsignc2waves,
    ("sig", "nc2diwasp"): runsignc2diwasp,
    ("vec", "dat2cdf"): runvecdat2cdf,
    ("vec", "cdf2nc"): runveccdf2nc,
    ("vec", "nc2waves"): runvecnc2waves,
    ("eco", "csv2cdf"): runecocsv2cdf,
    ("eco", "cdf2nc"): runecocdf2nc,
    ("eofe", "log2cdf"): runeofelog2cdf,
    ("eofe", "cdf2nc"): runeofecdf2nc,
    ("exo", "csv2cdf"): runexocsv2cdf,
    ("exo", "cdf2nc"): runexocdf2nc,
    ("glx", "dat2cdf"): runglxdat2cdf,
    ("glx", "cdf2nc"): runglxcdf2nc,
    ("glx", "nc2waves"): runglxnc2waves,
    ("hobo", "csv2cdf"): runhobocsv2cdf,
    ("hobo", "cdf2nc"): runhobocdf2nc,
    ("iq", "mat2cdf"): runiqmat2cdf,
    ("iq", "cdf2nc"): runiqcdf2nc,
    ("lisst", "csv2cdf"): runlisstcsv2cdf,
    ("lisst", "cdf2nc"): runlisstcdf2nc,
    ("mc", "asc2cdf"): runmcasc2cdf,
    ("mc", "cdf2nc"): runmccdf2nc,
    ("rdi", "mat2cdf"): lambda args: stglib.rdi.mat2cdf.mat_to_cdf(get_metadata(args)),
    ("rdi", "cdf2nc"): runrdicdf2nc,
    ("sgtid", "tid2cdf"): lambda args: stglib.sg.tid2cdf.tid_to_cdf(get_metadata(args)),
    ("sgtid", "cdf2nc"): lambda args: run_cdf_to_nc(stglib.sg.cdf2nc.cdf_to_nc, args),
    ("sgwvs", "wb2cdf"): lambda args: stglib.sg.wvswb2cdf.wb_to_cdf(get_metadata(args)),
    ("sgwvs", "cdf2nc"): lambda args: run_cdf_to_nc(
        stglib.sg.wvscdf2nc.cdf_to_nc, args
    ),
    ("sgwvs", "nc2waves"): lambda args: stglib.sg.wvsnc2waves.nc_to_waves(args.ncname),
    ("tb", "csv2cdf"): lambda args: stglib.tb.txt_to_cdf(get_metadata(args)),
    ("tb", "cdf2nc"): lambda args: run_cdf_to_nc(stglib.tb.cdf_to_nc, args),
    # TruBlue uses RSK nc2waves
    ("tb", "nc2waves"): lambda args: stglib.rsk.nc2waves.nc_to_waves(args.ncname),
    ("tcm", "csv2cdf"): runtcmcsv2cdf,
    ("tcm", "cdf2nc"): runtcmcdf2nc,
    ("troll", "csv2cdf"): runtrollcsv2cdf,
    ("troll", "cdf2nc"): runtrollcdf2nc,
    ("met", "csv2cdf"): runmetcsv2cdf,
    ("met", "cdf2nc"): runmetcdf2nc,
    ("son", "raw2cdf"): lambda args: stglib.son.raw2cdf.file81R_to_cdf(
        get_metadata(args)
    ),
    ("son", "cdf2nc"): lambda args: run_cdf_to_nc(stglib.son.cdf2nc.cdf_to_nc, args),
    ("son", "nc2xy"): lambda args: stglib.son.nc2xy.nc_to_xy(args.ncname),
    ("mar", "csv2cdf"): lambda args: stglib.mar.csv_to_cdf(get_metadata(args)),
    ("mar", "cdf2nc"): lambda args: run_cdf_to_nc(stglib.mar.cdf_to_nc, args),
}


def runots():
    args = stglib.cmd.runots_parser().parse_args()

    print(f"stglib {stglib.__version__}")

    key = (args.instrument, args.step)
    fn = _REGISTRY.get(key)
    if fn is None:
        raise ValueError(
            f"Unknown instrument/step combination: instrument={args.instrument!r}, step={args.step!r}"
        )
    fn(args)
