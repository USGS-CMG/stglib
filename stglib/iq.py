from __future__ import division, print_function
import scipy.io as spio
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import numpy as np
from . import core

def read_iq(filnam, start, stop, freq):
    iqmat = core.utils.loadmat(filnam)

    time = pd.date_range(start, stop, freq=freq)
    ds = {}
    ds['time'] = xr.DataArray(time, dims='time')
    attrs = {}
    for k in iqmat:
        if '__' not in k:
            if len(np.ravel(iqmat[k])) == len(time):
                ds[k] = xr.DataArray(np.ravel(iqmat[k]), dims='time')
                if k in iqmat['Data_Units']:
                    ds[k].attrs['units'] = iqmat['Data_Units'][k]
            elif '_2_' in k or '_3_' in k:
                ds[k] = xr.DataArray(iqmat[k], dims=('time', 'bin_beam_2_3'))
            elif '_0_' in k or '_1_' in k:
                ds[k] = xr.DataArray(iqmat[k], dims=('time', 'bin_beam_0_1'))
            elif 'FlowData_Vel' in k or 'FlowData_SNR' in k:
                ds[k] = xr.DataArray(iqmat[k], dims=('time', 'velbeam'))
            elif 'FlowData_NoiseLevel' in k:
                ds[k] = xr.DataArray(iqmat[k], dims=('time', 'beam'))

    ds = xr.Dataset(ds)
    for k in iqmat['System_IqSetup']['basicSetup']:
        if 'spare' not in k:
            ds.attrs[k] = iqmat['System_IqSetup']['basicSetup'][k]
    for k in iqmat['System_Id']:
        ds.attrs[k] = iqmat['System_Id'][k]
    for k in iqmat['System_IqState']:
        if 'spare' not in k:
            ds.attrs[k] = iqmat['System_IqState'][k]

    return ds

def make_iq_plots(iq, directory='', savefig=True):
    plt.figure(figsize=(11,8.5))

    for n, var in enumerate(['FlowData_Depth', 'FlowData_Vel_Mean', 'FlowData_Flow'], start=1):
        plt.subplot(3,1,n)
        plt.plot(iq['time'], iq[var])
        plt.ylabel(var + ' [' + iq[var].attrs['units'] + ']')
        # if n == 1:
            # plt.title(iq.attrs[''])

    if savefig:
        plt.savefig(directory + '/iq_stage_vel_flow.pdf')
    plt.show()


    # plt.figure(figsize=(11,8.5))
