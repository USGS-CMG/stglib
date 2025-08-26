import argparse

import stglib


def yamlarg(parser):
    parser.add_argument(
        "config", help="path to instrument configuration file (YAML formatted)"
    )


def gattsarg(parser):
    parser.add_argument(
        "gatts", help="path to global attributes file (gatts formatted)"
    )


def cdfarg(parser):
    parser.add_argument("cdfname", help="raw .cdf filename")


def atmarg(parser):
    parser.add_argument(
        "--atmpres", help="path to cdf file containing atmospheric pressure data"
    )


def hgtarg(parser):
    parser.add_argument(
        "--height", help="path to nc file containing height above seabed data"
    )


def swtarg(parser):
    parser.add_argument(
        "--salwtemp",
        help="path to nc file containing salinity and water temperature data",
    )


def ncarg(parser):
    parser.add_argument("ncname", help="processed .nc filename")


def addcdf2nc(instsp, description="Convert raw .cdf to clean .nc"):
    instsp.add_parser(
        "cdf2nc", parents=[cdf2nc_parser()], add_help=False, description=description
    )


def addnc2waves(instsp):
    instsp.add_parser("nc2waves", parents=[nc2waves_parser()], add_help=False)


def addnc2xy(instsp):
    instsp.add_parser("nc2xy", parents=[nc2xy_parser()], add_help=False)


def addnc2diwasp(instsp):
    instsp.add_parser("nc2diwasp", parents=[nc2diwasp_parser()], add_help=False)


def addinst2cdf(instsp, action, description="Convert instrument data to raw .cdf"):
    instsp.add_parser(
        action, parents=[inst2cdf_parser()], add_help=False, description=description
    )


def add_instrument(subparsers, instrument, description=None):
    inst = subparsers.add_parser(instrument, description=description)
    instsp = inst.add_subparsers(
        title="Steps",
        required=True,
        dest="step",
        description="Specify one of the steps in the list below",
    )
    return instsp


def runots_parser():
    parser = argparse.ArgumentParser(
        description="Run USGS CMHRP ocean time-series data processing system."
    )

    parser.add_argument(
        "-v", "--version", action="version", version=stglib.utils.get_detailed_version()
    )

    subparsers = parser.add_subparsers(
        title="Instruments",
        required=True,
        dest="instrument",
        description="Specify one of the instruments in the list below",
    )

    instsp = add_instrument(subparsers, "abss", "AQUAscat1000R")
    addinst2cdf(instsp, "mat2cdf")
    addcdf2nc(instsp)

    instsp = add_instrument(subparsers, "aqd", "Aquadopp (currents)")
    addinst2cdf(instsp, "hdr2cdf")
    addcdf2nc(instsp)

    instsp = add_instrument(subparsers, "aqdhr", "Aquadopp HR")
    addinst2cdf(instsp, "hdr2cdf")
    addcdf2nc(instsp)

    instsp = add_instrument(subparsers, "aqdwvs", "Aquadopp Waves")
    addinst2cdf(instsp, "wad2cdf")
    addcdf2nc(instsp)
    addnc2waves(instsp)

    instsp = add_instrument(subparsers, "eco", "WET Labs ECO")
    addinst2cdf(instsp, "csv2cdf")
    addcdf2nc(instsp)

    instsp = add_instrument(subparsers, "eofe", "EofE ECHOLOGGER")
    addinst2cdf(instsp, "log2cdf")
    addcdf2nc(instsp)

    instsp = add_instrument(subparsers, "exo", "YSI EXO")
    addinst2cdf(instsp, "csv2cdf")
    addcdf2nc(instsp)

    instsp = add_instrument(subparsers, "glx", "Geolux Wave Radar")
    addinst2cdf(instsp, "dat2cdf")
    addcdf2nc(instsp)
    addnc2waves(instsp)

    instsp = add_instrument(subparsers, "hobo", "Onset HOBO")
    addinst2cdf(instsp, "csv2cdf")
    addcdf2nc(instsp)

    instsp = add_instrument(subparsers, "iq", "SonTek IQ")
    addinst2cdf(instsp, "mat2cdf")
    addcdf2nc(instsp)

    instsp = add_instrument(subparsers, "lisst", "Sequoia Scientific LISST")
    addinst2cdf(instsp, "csv2cdf")
    addcdf2nc(instsp)

    instsp = add_instrument(subparsers, "mc", "Seabird MicroCAT")
    addinst2cdf(instsp, "asc2cdf")
    addcdf2nc(instsp)

    instsp = add_instrument(subparsers, "rbr", "RBR")
    addinst2cdf(instsp, "csv2cdf")
    addcdf2nc(instsp)
    addnc2waves(instsp)
    addnc2diwasp(instsp)

    instsp = add_instrument(subparsers, "rdi", "Teledyne RDI ADCP")
    addinst2cdf(instsp, "mat2cdf")
    addcdf2nc(instsp)

    instsp = add_instrument(subparsers, "rsk", "RBR")
    addinst2cdf(instsp, "csv2cdf")
    addcdf2nc(instsp)
    addnc2waves(instsp)
    addnc2diwasp(instsp)

    instsp = add_instrument(subparsers, "sgtid", "Seabird Seagauge Tides")
    addinst2cdf(instsp, "tid2cdf")
    addcdf2nc(instsp)

    instsp = add_instrument(subparsers, "sgwvs", "Seabird Seagauge Waves")
    addinst2cdf(instsp, "wb2cdf")
    addcdf2nc(instsp)
    addnc2waves(instsp)

    instsp = add_instrument(subparsers, "sig", "Nortek Signature")
    addinst2cdf(instsp, "mat2cdf")
    addcdf2nc(instsp)
    addnc2waves(instsp)
    addnc2diwasp(instsp)

    instsp = add_instrument(subparsers, "tb", "TruBlue")
    addinst2cdf(instsp, "csv2cdf")
    addcdf2nc(instsp)
    addnc2waves(instsp)

    instsp = add_instrument(subparsers, "tcm", "Lowell Tilt Current Meter")
    addinst2cdf(instsp, "csv2cdf")
    addcdf2nc(instsp)

    instsp = add_instrument(subparsers, "troll", "AquaTROLL")
    addinst2cdf(instsp, "csv2cdf")
    addcdf2nc(instsp)

    instsp = add_instrument(subparsers, "vec", "Nortek Vector")
    addinst2cdf(instsp, "dat2cdf")
    addcdf2nc(instsp)
    addnc2waves(instsp)

    instsp = add_instrument(subparsers, "wvs", "Aquadopp Waves")
    addinst2cdf(instsp, "wad2cdf")
    addcdf2nc(instsp)
    addnc2waves(instsp)

    instsp = add_instrument(subparsers, "wxt", "Vaisala WXT")
    addinst2cdf(instsp, "csv2cdf")
    addcdf2nc(instsp)

    instsp = add_instrument(subparsers, "son", "Imagenex sonar")
    addinst2cdf(instsp, "raw2cdf")
    addcdf2nc(instsp)
    addnc2xy(instsp)

    return parser


def inst2cdf_parser(description="Convert instrument files to raw .cdf format"):
    """generic parser for instrument data to raw .cdf; requires gatts and yaml"""
    # description = "Convert Aquadopp text files to raw .cdf format. Run this script from the directory containing Aquadopp files."
    parser = argparse.ArgumentParser(description=description)
    gattsarg(parser)
    yamlarg(parser)

    return parser


def cdf2nc_parser(
    description="Convert raw .cdf format to processed .nc files, optionally compensating for atmospheric pressure",
    atmpres=True,
    height=True,
    salwtemp=True,
):
    """generic parser for raw .cdf format to processed .nc files, optionally compensating for atmospheric pressure"""
    parser = argparse.ArgumentParser(description=description)
    cdfarg(parser)
    if atmpres:
        atmarg(parser)
    if height:
        hgtarg(parser)
    if salwtemp:
        swtarg(parser)

    return parser


def nc2waves_parser(description="Generate wave-statistics file", salwtemp=True):
    """generic parser for processed .nc to wave statistics"""
    parser = argparse.ArgumentParser(description=description)
    ncarg(parser)
    if salwtemp:
        swtarg(parser)

    return parser


def nc2xy_parser(description="Convert polar to cartesian coordinates"):
    """generic parser for processing sonar .nc to xy coordinates"""
    parser = argparse.ArgumentParser(description=description)
    ncarg(parser)

    return parser


def nc2diwasp_parser(description="Generate DIWASP wave-statistics file", salwtemp=True):
    """generic parser for processed .nc to wave statistics"""
    parser = argparse.ArgumentParser(description=description)
    ncarg(parser)
    if salwtemp:
        swtarg(parser)

    return parser


def aqdturnaround_parser():
    description = "Create Aquadopp turnaround plots. Run this script from the directory containing Aquadopp files."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("basefile", help="basename of AQD file (without extension)")
    parser.add_argument(
        "--orientation",
        default="UP",
        help="instrument orientation (UP/DOWN). Default UP",
    )

    return parser


def exoturnaround_parser():
    description = "Create YSI EXO turnaround plots. Run this script from the directory containing EXO files."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "basefile", help="basename of EXO .csv file (without extension)"
    )
    parser.add_argument(
        "--skiprows",
        default=25,
        type=int,
        help="Number of header rows to skip. Default 25",
    )

    return parser


def sigcdf2nc_parser():
    description = "Convert raw Signature .cdf format to processed .nc files, optionally compensating for atmospheric pressure"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("cdfname", help="raw .CDF filename(s)", nargs="+")
    atmarg(parser)

    return parser
