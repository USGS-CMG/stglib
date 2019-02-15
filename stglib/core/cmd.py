from __future__ import division, print_function
import argparse


def yamlarg(parser):
    parser.add_argument('config', help=('path to instrument configuration file'
                                        ' (YAML formatted)'))


def gattsarg(parser):
    parser.add_argument('gatts', help=('path to global attributes file (gatts '
                                       'formatted)'))


def aqdhdr2cdf_parser():
    description = ('Convert Aquadopp text files to raw .cdf format. Run this '
                   'script from the directory containing Aquadopp files')
    parser = argparse.ArgumentParser(description=description)
    gattsarg(parser)
    yamlarg(parser)

    return parser


def aqdcdf2nc_parser():
    description = ('Convert raw Aquadopp .cdf format to processed .nc files, '
                   'optionally compensating for atmospheric pressure')
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('cdfname', help='raw .CDF filename')
    parser.add_argument('--atmpres', help=('path to cdf file containing '
                                           'atmopsheric pressure data'))

    return parser


def wvswad2cdf_parser():
    description = ('Convert Aquadopp .wad wave files to raw .cdf format. '
                   'Run this script from the directory containing Aquadopp '
                   'files')
    parser = argparse.ArgumentParser(description=description)
    gattsarg(parser)
    yamlarg(parser)

    return parser


def wvscdf2nc_parser():
    description = 'Convert raw Aquadopp .cdf wave files to processed .nc files'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('cdfname', help='raw .CDF filename')
    parser.add_argument('--atmpres', help=('path to cdf file containing '
                                           'atmopsheric pressure data'))

    return parser


def wvsnc2diwasp_parser():
    description = 'Convert processed Aquadopp waves .nc files using DIWASP'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('ncname', help='processed .nc filename')

    return parser


def wvsnc2waves_parser():
    description = 'Generate waves statistics file'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('ncname', help='processed .nc filename')

    return parser

# MM 2/12/2019 removed mention of DWave so any RBR instrument will applys
def rskrsk2cdf_parser():
    description = ('Convert raw RBR files (.rsk) to raw .cdf format. '
                   'Run this script from the directory '
                   'containing .rsk files')
    parser = argparse.ArgumentParser(description=description)
    gattsarg(parser)
    yamlarg(parser)

    return parser


def rskcdf2nc_parser():
    description = ('Convert raw RBR .cdf format to processed .nc files,'
                   ' optionally compensating for atmospheric pressure')
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('cdfname', help='raw .CDF filename')
    parser.add_argument('--atmpres', help=('path to cdf file containing '
                                           'atmopsheric pressure data'))

    return parser


def rsknc2diwasp_parser():
    description = 'Convert processed .nc files using DIWASP'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('ncname', help='processed .nc filename')

    return parser


def rsknc2waves_parser():
    description = 'Generate waves statistics file'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('ncname', help='processed .nc filename')

    return parser


def hlmcsv2cdf_parser():
    description = ('Convert HOBO pressure sensor .csv file to raw .cdf format.'
                   ' Run this script from the directory containing IQ file')
    parser = argparse.ArgumentParser(description=description)
    gattsarg(parser)
    yamlarg(parser)

    return parser


def hlmcdf2nc_parser():
    description = 'Convert raw HOBO .cdf format to processed .nc files'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('cdfname', help='raw .CDF filename')

    return parser


def exocsv2cdf_parser():
    description = ('Convert EXO .csv file to raw .cdf format. Run this script '
                   'from the directory containing EXO file')
    parser = argparse.ArgumentParser(description=description)
    gattsarg(parser)
    yamlarg(parser)

    return parser


def exocdf2nc_parser():
    description = 'Convert raw EXO .cdf format to processed .nc files'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('cdfname', help='raw .CDF filename')
    parser.add_argument('--atmpres', help=('path to cdf file containing '
                                           'atmopsheric pressure data'))

    return parser


def ecolog2cdf_parser():
    description = ('Convert WET Labs ECO file to raw .cdf format. Run this '
                   'script from the directory containing ECO file')
    parser = argparse.ArgumentParser(description=description)
    gattsarg(parser)
    yamlarg(parser)

    return parser


def ecocdf2nc_parser():
    description = 'Convert raw ECO .cdf format to processed .nc files'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('cdfname', help='raw .CDF filename')

    return parser


def aqdturnaround_parser():
    description = ('Create Aquadopp turnaround plots. Run this script from '
                   'the directory containing Aquadopp files')
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('basefile',
                        help='basename of AQD file (without extension)')
    parser.add_argument('--orientation',
                        default='UP',
                        help='instrument orientation (UP/DOWN). Default UP')

    return parser


def exoturnaround_parser():
    description = ('Create YSI EXO turnaround plots. Run this script from '
                   'the directory containing EXO files')
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('basefile',
                        help='basename of EXO .csv file (without extension)')
    parser.add_argument('--skiprows',
                        default=25,
                        type=int,
                        help='Number of header rows to skip. Default 25')

    return parser
