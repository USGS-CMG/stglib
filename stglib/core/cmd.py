from __future__ import division, print_function
import argparse

def aqdhdr2cdf_parser():
    parser = argparse.ArgumentParser(description='Convert Aquadopp text files to raw .cdf format. Run this script from the directory containing Aquadopp files')
    parser.add_argument('gatts', help='path to global attributes file (gatts formatted)')
    parser.add_argument('config', help='path to ancillary config file (YAML formatted)')

    return parser

def aqdcdf2nc_parser():
    parser = argparse.ArgumentParser(description='Convert raw Aquadopp .cdf format to processed .nc files')
    parser.add_argument('cdfname', help='raw .CDF filename')
    parser.add_argument('--atmpres', help='path to cdf file containing atmopsheric pressure data')

    return parser
