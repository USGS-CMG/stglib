from __future__ import division, print_function
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
from . import core


def read_iq(filnam):
    """Read SonTek IQ data which has been exported as a Matlab .mat file from IQ
    software into an xarray Dataset

    Parameters
    ----------
    filnam : string
        The SonTek .mat filename

    Returns
    -------
    xarray.Dataset
        An xarray Dataset of the IQ data
    """

    iqmat = core.utils.loadmat(filnam)
    # offset = iqmat['FlowSubData_PrfHeader_0_BlankingDistance']
    # beamdist_0 = np.linspace(offset, offset + \
    # 100*iqmat['FlowSubData_PrfHeader_0_CellSize'], 100)
    ds = {}

    ds['time'] = xr.DataArray(
        iqmat['FlowData_SampleTime'],
        attrs={'standard_name': 'time',
               'axis': 'T',
               # per email from SonTek
               'units': 'microseconds since 2000-01-01 00:00:00',
               'calendar': 'proleptic_gregorian'}, dims='time')

    ds['velbeam'] = xr.DataArray([1, 2, 3, 4], dims='velbeam')
    ds['beam'] = xr.DataArray([1, 2, 3, 4, 5], dims='beam')
    # ds['beamdist_0'] = xr.DataArray(beamdist_0, dims='beamdist_0')
    # attrs = {}

    # need to do this because sometimes the flowsubdata and profile data is
    # one burst longer
    timelen = len(ds['time'])

    for k in iqmat:
        if '__' not in k and 'FlowSubData' not in k:
            print(k, np.shape(iqmat[k]))
            if len(np.ravel(iqmat[k])) == len(ds['time']):
                ds[k] = xr.DataArray(np.ravel(iqmat[k]), dims='time')
                if k in iqmat['Data_Units']:
                    ds[k].attrs['units'] = iqmat['Data_Units'][k]
            elif '_2_' in k or '_3_' in k:
                ds[k] = xr.DataArray(iqmat[k][0:timelen,:],
                                     dims=('time', 'beamdist_2_3'))
            elif '_0_' in k or '_1_' in k:
                ds[k] = xr.DataArray(iqmat[k][0:timelen,:],
                                     dims=('time', 'beamdist_0_1'))
            elif 'FlowData_Vel' in k or 'FlowData_SNR' in k:
                ds[k] = xr.DataArray(iqmat[k][0:timelen,:],
                                     dims=('time', 'velbeam'))
            elif 'FlowData_NoiseLevel' in k:
                ds[k] = xr.DataArray(iqmat[k][0:timelen,:],
                                     dims=('time', 'beam'))

    ds = xr.Dataset(ds)
    for k in iqmat['System_IqSetup']['basicSetup']:
        if 'spare' not in k:
            ds.attrs[k] = iqmat['System_IqSetup']['basicSetup'][k]
    for k in iqmat['System_Id']:
        ds.attrs[k] = iqmat['System_Id'][k]
    for k in iqmat['System_IqState']:
        if 'spare' not in k:
            ds.attrs[k] = iqmat['System_IqState'][k]

    return xr.decode_cf(ds)


def clean_iq(iq):
    """
    Preliminary data cleaning when SNR < 0
    """

    iq['FlowData_Vel_Mean'].values[iq['FlowData_Vel_Mean'] < -214748] = np.nan
    iq['FlowData_Vel'].values[iq['FlowData_Vel'] == -214748368] = np.nan
    for bm in range(4):
        pr = 'Profile_' + str(bm) + '_Vel'
        iq[pr].values[iq[pr] == -214748368] = np.nan

    return iq


def vel_to_ms(iq):
    """
    Convert velocity data from mm/s to m/s
    """

    for var in ['FlowData_Vel_Mean', 'FlowData_Vel']:
        iq[var] = iq[var] / 1000

    return iq


def make_beamdist(iq):
    """
    Generate physical coordinates to pair with the logical beamdist coordinates
    """
    for bm in range(4):
        if bm < 2:
            bdname = 'beamdist_0_1'
        else:
            bdname = 'beamdist_2_3'

        r = range(len(iq[bdname]))

        time = np.tile(iq['time'], (len(iq[bdname]), 1)).transpose()

        cells = np.zeros(np.shape(iq['Profile_' + str(bm) + '_Vel']))
        fsdbd = 'FlowSubData_PrfHeader_' + str(bm) + '_BlankingDistance'
        fsdcs = 'FlowSubData_PrfHeader_' + str(bm) + '_CellSize'
        for n in range(len(iq['time'])):
            cells[n, :] = iq[fsdbd][n].values + r * iq[fsdcs][n].values
        iq['cells_' + str(bm)] = xr.DataArray(cells, dims=('time', bdname))
        iq['time_' + str(bm)] = xr.DataArray(time, dims=('time', bdname))
        iq.set_coords(['cells_' + str(bm),
                       'time_' + str(bm)],
                       inplace=True)

    return iq


def make_iq_plots(iq, directory='', savefig=False):
    """
    Make IQ turnaround plots
    """

    plt.figure(figsize=(11, 8.5))

    for n, var in enumerate(['FlowData_Depth',
                             'FlowData_Vel_Mean',
                             'FlowData_Flow'], start=1):
        plt.subplot(3, 1, n)
        plt.plot(iq['time'], iq[var])
        plt.ylabel(var + ' [' + iq[var].attrs['units'] + ']')

    if savefig:
        plt.savefig(directory + '/iq_stage_vel_flow.pdf')
    plt.show()
