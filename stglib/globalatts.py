from __future__ import division, print_function
import csv

def read_globalatts(fname):
    """
    Read global attributes file (glob_attxxxx.txt) and create metadata structure
    """

    metadata = {}

    with open(fname, 'r') as csvfile:
        a = csv.reader(csvfile, delimiter=';')

        for row in a:
            if row[0] == 'MOORING':
                metadata[row[0].strip()] = row[1].strip()
            else:
                metadata[row[0].strip()] = str2num(row[1].strip())

        return metadata

def str2num(s):
    """
    Convert string to float if possible
    """

    try:
        float(s)
        return float(s)
    except ValueError:
        return s
