from __future__ import division, print_function
import argparse

def yamlarg(parser):
    parser.add_argument('config', help='path to instrument configuration file (YAML formatted)')

def gattsarg(parser):
    parser.add_argument('gatts', help='path to global attributes file (gatts formatted)')

def aqdhdr2cdf_parser():
    parser = argparse.ArgumentParser(description='Convert Aquadopp text files to raw .cdf format. Run this script from the directory containing Aquadopp files')
    gattsarg(parser)
    yamlarg(parser)

    return parser

def aqdhdr2cdf_parse_args():
    parser = aqdhdr2cdf_parser()

    return parser.parse_args()

def aqdcdf2nc_parser():
    parser = argparse.ArgumentParser(description='Convert raw Aquadopp .cdf format to processed .nc files, optionally compensating for atmospheric pressure')
    parser.add_argument('cdfname', help='raw .CDF filename')
    parser.add_argument('--atmpres', help='path to cdf file containing atmopsheric pressure data')

    return parser

def aqdcdf2nc_parse_args():
    parser = aqdcdf2nc_parser()

    return parser.parse_args()

def rskrsk2cdf_parser():
    parser = argparse.ArgumentParser(description='Convert raw RBR d|wave files (.rsk) to raw .cdf format. Run this script from the directory containing d|wave files')
    gattsarg(parser)
    yamlarg(parser)

    return parser

def rskrsk2cdf_parse_args():
    parser = rskrsk2cdf_parser()

    return parser.parse_args()

def rskcdf2nc_parser():
    parser = argparse.ArgumentParser(description='Convert raw RBR d|wave .cdf format to processed .nc files, optionally compensating for atmospheric pressure')
    parser.add_argument('cdfname', help='raw .CDF filename')
    parser.add_argument('--atmpres', help='path to cdf file containing atmopsheric pressure data')

    return parser

def rskcdf2nc_parse_args():
    parser = rskcdf2nc_parser()

    return parser.parse_args()

def rsknc2diwasp_parser():
    parser = argparse.ArgumentParser(description='Convert processed .nc files using DIWASP')
    parser.add_argument('ncname', help='processed .nc filename')

    return parser

def rsknc2diwasp_parse_args():
    parser = rsknc2diwasp_parser()

    return parser.parse_args()
