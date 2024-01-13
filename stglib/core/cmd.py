import argparse


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
        "--atmpres", help="path to cdf file containing atmopsheric pressure data"
    )


def ncarg(parser):
    parser.add_argument("ncname", help="processed .nc filename")


def aqdhdr2cdf_parser():
    description = "Convert Aquadopp text files to raw .cdf format. Run this script from the directory containing Aquadopp files."
    parser = argparse.ArgumentParser(description=description)
    gattsarg(parser)
    yamlarg(parser)

    return parser


def aqdcdf2nc_parser():
    description = "Convert raw Aquadopp .cdf format to processed .nc files, optionally compensating for atmospheric pressure"
    parser = argparse.ArgumentParser(description=description)
    cdfarg(parser)
    atmarg(parser)

    return parser


def wvswad2cdf_parser():
    description = "Convert Aquadopp .wad wave files to raw .cdf format. Run this script from the directory containing Aquadopp files."
    parser = argparse.ArgumentParser(description=description)
    gattsarg(parser)
    yamlarg(parser)

    return parser


def wvscdf2nc_parser():
    description = "Convert raw Aquadopp .cdf wave files to processed .nc files"
    parser = argparse.ArgumentParser(description=description)
    cdfarg(parser)
    atmarg(parser)

    return parser


def wvsnc2diwasp_parser():
    description = "Convert processed Aquadopp waves .nc files using DIWASP"
    parser = argparse.ArgumentParser(description=description)
    ncarg(parser)

    return parser


def wvsnc2waves_parser():
    description = "Generate waves statistics file"
    parser = argparse.ArgumentParser(description=description)
    ncarg(parser)

    return parser


def rskcsv2cdf_parser():
    description = "Convert exported RBR csv files to raw .cdf format. Run this script from the directory containing the files."
    parser = argparse.ArgumentParser(description=description)
    gattsarg(parser)
    yamlarg(parser)

    return parser


def rskrsk2cdf_parser():
    description = "Convert raw RBR files (.rsk) to raw .cdf format. Run this script from the directory containing the files."
    parser = argparse.ArgumentParser(description=description)
    gattsarg(parser)
    yamlarg(parser)

    return parser


def rskcdf2nc_parser():
    description = "Convert raw RBR d|wave .cdf format to processed .nc files, optionally compensating for atmospheric pressure"
    parser = argparse.ArgumentParser(description=description)
    cdfarg(parser)
    atmarg(parser)

    return parser


def rsknc2diwasp_parser():
    description = "Convert processed .nc files using DIWASP"
    parser = argparse.ArgumentParser(description=description)
    ncarg(parser)

    return parser


def rsknc2waves_parser():
    description = "Generate waves statistics file"
    parser = argparse.ArgumentParser(description=description)
    ncarg(parser)

    return parser


def hobocsv2cdf_parser():
    description = "Convert HOBO pressure sensor .csv file to raw .cdf format. Run this script from the directory containing HOBO file."
    parser = argparse.ArgumentParser(description=description)
    gattsarg(parser)
    yamlarg(parser)

    return parser


def hobocdf2nc_parser():
    description = "Convert raw HOBO .cdf format to processed .nc files"
    parser = argparse.ArgumentParser(description=description)
    cdfarg(parser)

    return parser


def iqmat2cdf_parser():
    description = "Convert SonTek IQ .mat file to raw .cdf format. Run this script from the directory containing IQ file."
    parser = argparse.ArgumentParser(description=description)
    gattsarg(parser)
    yamlarg(parser)

    return parser


def iqcdf2nc_parser():
    description = "Convert raw SonTek IQ .cdf format to processed .nc files"
    parser = argparse.ArgumentParser(description=description)
    cdfarg(parser)

    return parser


def exocsv2cdf_parser():
    description = "Convert EXO .csv file to raw .cdf format. Run this script from the directory containing EXO file."
    parser = argparse.ArgumentParser(description=description)
    gattsarg(parser)
    yamlarg(parser)

    return parser


def exocdf2nc_parser():
    description = "Convert raw EXO .cdf format to processed .nc files"
    parser = argparse.ArgumentParser(description=description)
    cdfarg(parser)
    atmarg(parser)

    return parser


def lisstcsv2cdf_parser():
    description = "Convert LISST .csv file to raw .cdf format. Run this script from the directory containing LISST file."
    parser = argparse.ArgumentParser(description=description)
    gattsarg(parser)
    yamlarg(parser)

    return parser


def lisstcdf2nc_parser():
    description = "Convert raw LISST .cdf format to processed .nc files"
    parser = argparse.ArgumentParser(description=description)
    cdfarg(parser)
    atmarg(parser)

    return parser


def trollcsv2cdf_parser():
    description = "Convert Aqua TROLL .csv file to raw .cdf format. Run this script from the directory containing .csv file."
    parser = argparse.ArgumentParser(description=description)
    gattsarg(parser)
    yamlarg(parser)

    return parser


def trollcdf2nc_parser():
    description = "Convert raw Aqua TROLL .cdf format to processed .nc files"
    parser = argparse.ArgumentParser(description=description)
    cdfarg(parser)

    return parser


def ecolog2cdf_parser():
    description = "Convert WET Labs ECO file to raw .cdf format. Run this script from the directory containing ECO file."
    parser = argparse.ArgumentParser(description=description)
    gattsarg(parser)
    yamlarg(parser)

    return parser


def ecocdf2nc_parser():
    description = "Convert raw ECO .cdf format to processed .nc files"
    parser = argparse.ArgumentParser(description=description)
    cdfarg(parser)

    return parser


def rdiraw2cdf_parser():
    description = "Convert RDI raw binary files to raw .cdf format. Run this script from the directory containing RDI files."
    parser = argparse.ArgumentParser(description=description)
    gattsarg(parser)
    yamlarg(parser)

    return parser


def rdicdf2nc_parser():
    description = "Convert raw RDI .cdf format to processed .nc files, optionally compensating for atmospheric pressure"
    parser = argparse.ArgumentParser(description=description)
    cdfarg(parser)
    atmarg(parser)

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


def wxtcsv2cdf_parser():
    description = "Convert Vaisala WXT met .csv file to raw .cdf format. Run this script from the directory containing Vaisala .csv file."
    parser = argparse.ArgumentParser(description=description)
    gattsarg(parser)
    yamlarg(parser)

    return parser


def wxtcdf2nc_parser():
    description = "Convert raw Vaisala WXT .cdf format to processed .nc files"
    parser = argparse.ArgumentParser(description=description)
    cdfarg(parser)

    return parser


def eofelog2cdf_parser():
    description = "Convert EofE echologger .log file to raw .cdf format. Run this script from the directory containing ea400 echologger .log file."
    parser = argparse.ArgumentParser(description=description)
    gattsarg(parser)
    yamlarg(parser)

    return parser


def eofecdf2nc_parser():
    description = "Convert raw echologger .cdf format to processed .nc files"
    parser = argparse.ArgumentParser(description=description)
    cdfarg(parser)

    return parser


def vechdr2cdf_parser():
    description = "Convert Vector text files to raw .cdf format. Run this script from the directory containing Vector files."
    parser = argparse.ArgumentParser(description=description)
    gattsarg(parser)
    yamlarg(parser)

    return parser


def veccdf2nc_parser():
    description = "Convert raw Vector .cdf format to processed .nc files, optionally compensating for atmospheric pressure"
    parser = argparse.ArgumentParser(description=description)
    cdfarg(parser)
    atmarg(parser)

    return parser


def sigmat2cdf_parser():
    description = "Convert Signature files exported in Matlab format to raw .cdf format. Run this script from the directory containing Signature files."
    parser = argparse.ArgumentParser(description=description)
    gattsarg(parser)
    yamlarg(parser)

    return parser


def sigraw2cdf_parser():
    description = "Convert raw Signature .ad2cp files to raw .cdf format. Run this script from the directory containing Signature files."
    parser = argparse.ArgumentParser(description=description)
    gattsarg(parser)
    yamlarg(parser)

    return parser


def sigcdf2nc_parser():
    description = "Convert raw Signature .cdf format to processed .nc files, optionally compensating for atmospheric pressure"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("cdfname", help="raw .CDF filename(s)", nargs="+")
    atmarg(parser)

    return parser


def sigdlfncdf2nc_parser():
    description = "Convert raw Signature .cdf format to processed .nc files, optionally compensating for atmospheric pressure"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("cdfname", help="raw .CDF filename(s)", nargs="+")
    atmarg(parser)

    return parser
