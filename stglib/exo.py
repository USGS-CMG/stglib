from __future__ import division, print_function
import pandas as pd
import xarray as xr
from .core import utils

def read_exo(filnam, skiprows=25):
    """Read data from a YSI EXO multiparameter sonde .xlsx file into an xarray
    Dataset.

    Parameters
    ----------
    filnam : string
        The filename
    skiprows : int, optional
        How many header rows to skip. Default 25

    Returns
    -------
    xarray.Dataset
        An xarray Dataset of the EXO data
    """

    exo = pd.read_excel(filnam,
                        skiprows=skiprows,
                        infer_datetime_format=True,
                        parse_dates=[['Date (MM/DD/YYYY)', 'Time (HH:MM:SS)']])
    exo.rename(columns={'Date (MM/DD/YYYY)_Time (HH:MM:SS)': 'time'}, inplace=True)
    exo.set_index('time', inplace=True)
    exo.rename(columns=lambda x: x.replace(' ', '_'), inplace=True)
    exo.rename(columns=lambda x: x.replace('/', '_per_'), inplace=True)
    exo['Press_dbar'] = exo['Press_psi_a'] * 0.689476

    return xr.Dataset(exo)

def xls_to_cdf(metadata):

    basefile = metadata['basefile']

    if 'skiprows' in metadata:
        ds = read_exo(basefile + '.xlsx', skiprows=metadata['skiprows'])
    else:
        ds = read_exo(basefile + '.xlsx')

    # write out metadata first, then deal exclusively with xarray attrs
    ds = utils.write_metadata(ds, metadata)

    del metadata

    ds = utils.create_epic_time(ds)

    # configure file
    cdf_filename = ds.attrs['filename'] + '-raw.cdf'

    ds.to_netcdf(cdf_filename, engine='netcdf4')

    print('Finished writing data to %s' % cdf_filename)

    return ds

def cdf_to_nc(cdf_filename, atmpres=False):
    """
    Load a "raw" .cdf file and generate a processed .nc file
    """

    # Load raw .cdf data
    ds = xr.open_dataset(cdf_filename, autoclose=True)

    # Clip data to in/out water times or via good_ens
    ds = utils.clip_ds(ds)

    # Rename some variables
    ds = ds.rename({'Press_dbar': 'P_1'})

    ds = ds.drop(['Press_psi_a'])

    if atmpres:
        print("Atmospherically correcting data")

        met = xr.open_dataset(atmpres, autoclose=True)
        # need to save attrs before the subtraction, otherwise they are lost
        attrs = ds['P_1'].attrs
        ds['P_1ac'] = ds['P_1'] - met['atmpres'] - met['atmpres'].offset
        print('Correcting using offset of %f' % met['atmpres'].offset)
        ds['P_1ac'].attrs = attrs
