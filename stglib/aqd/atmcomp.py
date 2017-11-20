from __future__ import division, print_function
import numpy as np
import calendar
import matplotlib.pyplot as plt

def tots(d):
    return calendar.timegm(d.timetuple())

def atmcomp(aqddatetime, aqdpress, metdatetime, metpress, offset=0):
    """
    Atmospheric pressure compensation of pressure records.
    Inputs:
        aqddatetime: array of datetimes from Aquadopp [decibars]
        aqdpress: array of pressure (depth) values from Aquadopp

        metdatetime: array of datetimes from met station
        metpress: array of pressure (atmospheric) *** NOTE [decibars] ***

        offset: offset for when sensor is out of water

    Outputs:
        aqdpress_ac: atmospherically corrected pressure record
    """

    ptime = np.array([tots(x) for x in aqddatetime])
    mtime = np.array([tots(x) for x in metdatetime])

    metinterp = np.interp(ptime, mtime, metpress)

    aqdpress_ac = aqdpress - (metinterp - offset)

    return aqdpress_ac

def plot_atmcomp(aqddatetime, aqdpress, aqdpress_ac, xlims=False, ylims=False):
    plt.figure(figsize=(12,8))
    plt.plot(aqddatetime, aqdpress)
    plt.plot(aqddatetime, aqdpress_ac)
    if xlims is not False:
        plt.xlim(xlims)
    if ylims is not False:
        plt.ylim(ylims)
    plt.show()
