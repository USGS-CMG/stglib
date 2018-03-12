from __future__ import division, print_function
import pandas as pd
import xarray as xr
from .core import utils

def read_exo(filnam, skiprows=25, encoding='utf-8'):
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

    exo = pd.read_csv(filnam,
                        skiprows=skiprows,
                        infer_datetime_format=True,
                        parse_dates=[['Date (MM/DD/YYYY)', 'Time (HH:MM:SS)']],
                        encoding=encoding)
    exo.rename(columns={'Date (MM/DD/YYYY)_Time (HH:MM:SS)': 'time'}, inplace=True)
    exo.set_index('time', inplace=True)
    exo.rename(columns=lambda x: x.replace(' ', '_'), inplace=True)
    exo.rename(columns=lambda x: x.replace('/', '_per_'), inplace=True)
    exo['Press_dbar'] = exo['Press_psi_a'] * 0.689476

    return xr.Dataset(exo)

def xls_to_cdf(metadata):

    basefile = metadata['basefile']

    try:
        ds = read_exo(basefile + '.csv', skiprows=metadata['skiprows'])
    except UnicodeDecodeError:
        # try reading as Mac OS Western for old versions of Mac Excel
        ds = read_exo(basefile + '.csv', skiprows=metadata['skiprows'], encoding='mac-roman')

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

    ds = ds_rename_vars(ds)

    # ds = ds_add_attrs(ds)

    ds = ds.drop(['Press_psi_a', 'Site_Name', 'Fault_Code', 'Time_(Fract._Sec)'])

    if atmpres:
        print("Atmospherically correcting data")

        met = xr.open_dataset(atmpres, autoclose=True)
        # need to save attrs before the subtraction, otherwise they are lost
        attrs = ds['P_1'].attrs
        ds['P_1ac'] = ds['P_1'] - met['atmpres'] - met['atmpres'].offset
        print('Correcting using offset of %f' % met['atmpres'].offset)
        ds['P_1ac'].attrs = attrs

    # assign min/max:
    ds = utils.add_min_max(ds)

    ds = utils.create_epic_time(ds)

    ds = ds_add_lat_lon(ds)

    ds = ds_add_attrs(ds)

    # Write to .nc file
    print("Writing cleaned/trimmed data to .nc file")
    nc_filename = ds.attrs['filename'] + '-a.nc'

    ds = utils.rename_time(ds)

    ds.to_netcdf(nc_filename)
    print('Done writing netCDF file', nc_filename)




def ds_rename_vars(ds):

    # set up dict of instrument -> EPIC variable names
    varnames = {'Press_dbar': 'P_1',
                    'Battery_V': 'Bat_106',
                    'fDOM_RFU': 'fDOMRFU',
                    'fDOM_QSU': 'fDOMQSU',
                    'Chlorophyll_RFU': 'CHLrfu', # capitalization based on Chincoteague names
                    'Chlorophyll_µg_per_L': 'Fch_906',
                    'BGA-PE_RFU': 'BGAPErfu',
                    'BGA-PE_µg_per_L': 'BGAPE',
                    'Temp_°C': 'T_28',
                    'Cond_mS_per_cm': 'C_51',
                    'SpCond_mS_per_cm': 'SpC_48',
                    'Sal_psu': 'S_41',
                    'ODO_%_sat': 'OST_62',
                    'ODO_mg_per_L': 'DO',
                    'Turbidity_NTU': 'Turb',
                    'pH': 'pH_159',
                    'pH_mV': 'pHmV',
                    'Depth_m': 'D_3'}

    # check to make sure they exist before trying to rename
    newvars = {}
    for k in varnames:
        if k in ds:
            newvars[k] = varnames[k]

    return ds.rename(newvars)

def ds_add_attrs(ds):
    # Update attributes for EPIC and STG compliance
    ds.lat.encoding['_FillValue'] = False
    ds.lon.encoding['_FillValue'] = False
    ds.time.encoding['_FillValue'] = False
    ds.epic_time.encoding['_FillValue'] = False
    ds.epic_time2.encoding['_FillValue'] = False

    ds['time'].attrs.update({'standard_name': 'time',
                             'axis': 'T'})

    ds['epic_time'].attrs.update({'units': 'True Julian Day',
                                  'type': 'EVEN',
                                  'epic_code': 624})

    ds['epic_time2'].attrs.update({'units': 'msec since 0:00 GMT',
                                   'type': 'EVEN',
                                   'epic_code': 624})

    ds['Bat_106'].attrs.update({'units': 'V',
                                'long_name': 'Battery voltage',
                                'epic_code': 106})

    ds['fDOMRFU'].attrs.update({'units': 'Relative fluorescence units (RFU)',
                                'long_name': 'Fluorescent dissolved organic matter'})

    ds['fDOMQSU'].attrs.update({'units': 'Quinine sulfate equivalent units (QSU)',
                                'long_name': 'Fluorescent dissolved organic matter'})

    ds['CHLrfu'].attrs.update({'units': 'Relative fluorescence units (RFU)',
                               'long_name': 'Chlorophyll A'})

    ds['Fch_906'].attrs.update({'units': 'ug/L',
                                'long_name': 'Chlorophyll A',
                                'epic_code': 906})

    ds['BGAPErfu'].attrs.update({'units': 'Relative fluorescence units (RFU)',
                                 'long_name': 'Blue green algae phycoerythrin'})

    ds['BGAPE'].attrs.update({'units': 'ug/L',
                              'long_name': 'Blue green algae phycoerythrin'})

    ds['T_28'].attrs.update({'units': 'C',
                            'long_name': 'Temperature',
                            'epic_code': 28})

    ds['C_51'].attrs.update({'units': 'S/m',
                            'long_name': 'Conductivity',
                            'epic_code': 51})
    ds['C_51'] = ds['C_51']/10 # convert from mS/cm to S/m

    ds['SpC_48'].attrs.update({'units': 'mS/cm',
                               'long_name': 'Conductivity',
                               'comment': 'Temperature compensated to 25 °C',
                               'epic_code': 48})

    ds['S_41'].attrs.update({'units': 'Practical salinity units (PSU)',
                             'long_name': 'Salinity',
                             'epic_code': 41})

    ds['OST_62'].attrs.update({'units': '%',
                              'long_name': 'Oxygen percent saturation',
                              'epic_code': 62})

    ds['DO'].attrs.update({'units': 'mg/L',
                           'long_name': 'Dissolved oxygen'})

    ds['Turb'].attrs.update({'units': 'Nephelometric turbidity units (NTU)',
                           'long_name': 'Turbidity'})

    ds['pH_159'].attrs.update({'units': '',
                               'long_name': 'pH',
                               'epic_code': 159})

    ds['P_1'].attrs.update({'units': 'mV',
                             'long_name': 'pH'})

    ds['P_1'].attrs.update({'units': 'dbar',
                            'long_name': 'Pressure',
                            'epic_code': 1})

    ds['D_3'].attrs.update({'units': 'm',
                            'long_name': 'Depth',
                            'epic_code': 3})

    if 'P_1ac' in ds:
        ds['P_1ac'].attrs.update({'units': 'dbar',
                                  'name': 'Pac',
                                  'long_name': 'Corrected pressure'})
        if 'P_1ac_note' in ds.attrs:
            ds['P_1ac'].attrs.update({'note': ds.attrs['P_1ac_note']})

    def add_attributes(var, dsattrs):
        var.attrs.update({
            'initial_instrument_height': dsattrs['initial_instrument_height'],
            # 'nominal_instrument_depth': dsattrs['nominal_instrument_depth'],
            'height_depth_units': 'm',
            })
        var.encoding['_FillValue'] = 1e35

    for var in ds.variables:
        if var not in ds.coords and var != 'time2':
            add_attributes(ds[var], ds.attrs)

    return ds

def ds_add_lat_lon(ds):
    ds['lat'] = xr.DataArray([ds.attrs['latitude']], dims=('lat'), name='lat')
    ds['lon'] = xr.DataArray([ds.attrs['longitude']], dims=('lon'), name='lon')

    return ds
