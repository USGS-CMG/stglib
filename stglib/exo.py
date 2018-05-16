from __future__ import division, print_function
import pandas as pd
import xarray as xr
import numpy as np
import scipy.signal
from .core import utils


def read_exo(filnam, skiprows=25, encoding='utf-8'):
    """Read data from a YSI EXO multiparameter sonde .csv file into an xarray
    Dataset.

    Parameters
    ----------
    filnam : string
        The filename
    skiprows : int, optional
        How many header rows to skip. Default 25
    encoding : string, optional
        File encoding. Default 'utf-8'

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
    exo.rename(columns={'Date (MM/DD/YYYY)_Time (HH:MM:SS)': 'time'},
               inplace=True)
    exo.set_index('time', inplace=True)
    exo.rename(columns=lambda x: x.replace(' ', '_'), inplace=True)
    exo.rename(columns=lambda x: x.replace('/', '_per_'), inplace=True)
    exo['Press_dbar'] = exo['Press_psi_a'] * 0.689476

    exo = xr.Dataset(exo)
    hdr = read_exo_header(filnam, encoding=encoding)
    exo.attrs['serial_number'] = hdr['serial_number']

    # Apply sensor serial numbers to each sensor
    for k in exo.variables:
        if 'fDOM' in k:
            exo[k].attrs['sensor_serial_number'] = (
                hdr['fDOM']['sensor_serial_number'])
        elif 'Chlorophyll' in k or 'BGA-PE' in k:
            exo[k].attrs['sensor_serial_number'] = (
                hdr['Total Algae BGA-PE']['sensor_serial_number'])
        elif 'Temp' in k or 'Cond' in k or 'Sal' in k:
            if 'Unknown CT' in hdr:
                exo[k].attrs['sensor_serial_number'] = (
                    hdr['Unknown CT']['sensor_serial_number'])
            elif 'Wiped CT' in hdr:
                exo[k].attrs['sensor_serial_number'] = (
                    hdr['Wiped CT']['sensor_serial_number'])
        elif 'ODO' in k:
            exo[k].attrs['sensor_serial_number'] = (
                hdr['Optical DO']['sensor_serial_number'])
        elif 'Turbidity' in k:
            exo[k].attrs['sensor_serial_number'] = (
                hdr['Turbidity']['sensor_serial_number'])
        elif 'pH' in k:
            exo[k].attrs['sensor_serial_number'] = (
                hdr['pH']['sensor_serial_number'])
        elif 'Press' in k or 'Depth' in k:
            exo[k].attrs['sensor_serial_number'] = (
                hdr['Depth Non-Vented 0-10m']['sensor_serial_number'])

    return exo


def csv_to_cdf(metadata):
    """
    Process EXO .csv file to a raw .cdf file
    """

    basefile = metadata['basefile']

    try:
        ds = read_exo(basefile + '.csv', skiprows=metadata['skiprows'])
    except UnicodeDecodeError:
        # try reading as Mac OS Western for old versions of Mac Excel
        ds = read_exo(basefile + '.csv',
                      skiprows=metadata['skiprows'],
                      encoding='mac-roman')

    # write out metadata first, then deal exclusively with xarray attrs
    ds = utils.write_metadata(ds, metadata)

    del metadata

    ds = utils.create_epic_time(ds)

    # configure file
    cdf_filename = ds.attrs['filename'] + '-raw.cdf'

    ds.to_netcdf(cdf_filename, unlimited_dims=['time'])

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

    ds = ds.drop(['Press_psi_a',
                  'Site_Name',
                  'Fault_Code',
                  'Time_(Fract._Sec)'])

    if atmpres:
        print("Atmospherically correcting data")

        met = xr.open_dataset(atmpres, autoclose=True)
        # need to save attrs before the subtraction, otherwise they are lost
        attrs = ds['P_1'].attrs
        ds['P_1ac'] = ds['P_1'] - met['atmpres'] - met['atmpres'].offset
        print('Correcting using offset of %f' % met['atmpres'].offset)
        ds['P_1ac'].attrs = attrs

    ds = exo_qaqc(ds)

    # assign min/max:
    ds = utils.add_min_max(ds)

    ds = utils.create_epic_time(ds)

    ds = ds_add_lat_lon(ds)

    ds = ds_add_attrs(ds)

    # Write to .nc file
    print("Writing cleaned/trimmed data to .nc file")
    nc_filename = ds.attrs['filename'] + '-a.nc'

    ds = utils.rename_time(ds)

    ds.to_netcdf(nc_filename, unlimited_dims=['time'])
    print('Done writing netCDF file', nc_filename)


def ds_rename_vars(ds):

    # set up dict of instrument -> EPIC variable names
    varnames = {'Press_dbar': 'P_1',
                'Battery_V': 'Bat_106',
                'fDOM_RFU': 'fDOMRFU',
                'fDOM_QSU': 'fDOMQSU',
                # capitalization based on Chincoteague names
                'Chlorophyll_RFU': 'CHLrfu',
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

    ds['fDOMRFU'].attrs.update({
        'units': 'Relative fluorescence units (RFU)',
        'long_name': 'Fluorescent dissolved organic matter'})

    ds['fDOMQSU'].attrs.update({
        'units': 'Quinine sulfate equivalent units (QSU)',
        'long_name': 'Fluorescent dissolved organic matter'})

    ds['CHLrfu'].attrs.update({
        'units': 'Relative fluorescence units (RFU)',
        'long_name': 'Chlorophyll A'})

    ds['Fch_906'].attrs.update({'units': 'ug/L',
                                'long_name': 'Chlorophyll A',
                                'epic_code': 906})

    ds['BGAPErfu'].attrs.update({
        'units': 'Relative fluorescence units (RFU)',
        'long_name': 'Blue green algae phycoerythrin'})

    ds['BGAPE'].attrs.update({
        'units': 'ug/L',
        'long_name': 'Blue green algae phycoerythrin'})

    ds['T_28'].attrs.update({'units': 'C',
                             'long_name': 'Temperature',
                             'epic_code': 28})

    ds['C_51'].attrs.update({'units': 'S/m',
                             'long_name': 'Conductivity',
                             'epic_code': 51})
    ds['C_51'] = ds['C_51']/10  # convert from mS/cm to S/m

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

    if 'pH_159' in ds.variables:
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


def read_exo_header(filnam, encoding='utf-8'):
    hdr = pd.read_csv(filnam, skiprows=None, encoding=encoding)
    hdr = pd.DataFrame(hdr.iloc[:, 0:4])
    # print(hdr)
    header = {}
    header['serial_number'] = (
        hdr[hdr['KOR Export File'] == 'Sonde ID'].values[0][1].split(' ')[1])
    for var in ['fDOM',
                'Total Algae BGA-PE',
                'Wiped CT',
                'Unknown CT',
                'Optical DO',
                'Turbidity',
                'pH',
                'Depth Non-Vented 0-10m']:
        vals = hdr[hdr['KOR Export File'] == var]
        if not vals.empty:
            header[var] = {}
            header[var]['sensor_serial_number'] = vals.values[0][1]
            header[var]['data_columns'] = (
                [int(x) for x in vals.values[0][3].split(';')])

    return header


def exo_qaqc(ds):
    """
    QA/QC
    Trim EXO data based on metadata
    """

    for var in ['C_51', 'SpC_48', 'S_41', 'Turb']:
        if var + '_min_diff' in ds.attrs:
            print('%s: Trimming using minimum diff of %f' %
                  (var, ds.attrs[var + '_min_diff']))
            ds[var][np.ediff1d(
                ds[var], to_begin=0) < ds.attrs[var + '_min_diff']] = np.nan

            notetxt = ('Values filled where data decreases by more than %f '
                       'units in a single time step. ' %
                       ds.attrs[var + '_min_diff'])

            if 'note' in ds[var].attrs:
                ds[var].attrs['note'] = notetxt + ds[var].attrs['note']
            else:
                ds[var].attrs.update({'note': notetxt})

        if var + '_max_diff' in ds.attrs:
            print('%s: Trimming using maximum diff of %f' %
                  (var, ds.attrs[var + '_max_diff']))
            ds[var][np.ediff1d(
                ds[var], to_begin=0) > ds.attrs[var + '_max_diff']] = np.nan

            notetxt = ('Values filled where data increases by more than %f '
                       'units in a single time step. ' %
                       ds.attrs[var + '_max_diff'])

            if 'note' in ds[var].attrs:
                ds[var].attrs['note'] = notetxt + ds[var].attrs['note']
            else:
                ds[var].attrs.update({'note': notetxt})

        if var + '_med_diff' in ds.attrs:
            if 'kernel_size' in ds.attrs:
                kernel_size = ds.attrs['kernel_size']
            else:
                kernel_size = 5
            print('%s: Trimming using %d-point median filter diff of %f' %
                  (var, kernel_size, ds.attrs[var + '_med_diff']))
            filtered = scipy.signal.medfilt(ds[var], kernel_size=kernel_size)
            bads = np.abs(ds[var] - filtered) > ds.attrs[var + '_med_diff']
            ds[var][bads] = np.nan

            notetxt = ('Values filled where difference between %d-point '
                       'median filter and original values is greater than '
                       '%f. ' %
                       (kernel_size, ds.attrs[var + '_med_diff']))

            if 'note' in ds[var].attrs:
                ds[var].attrs['note'] = notetxt + ds[var].attrs['note']
            else:
                ds[var].attrs.update({'note': notetxt})

        if var + '_med_diff_pct' in ds.attrs:
            if 'kernel_size' in ds.attrs:
                kernel_size = ds.attrs['kernel_size']
            else:
                kernel_size = 5
            print('%s: Trimming using %d-point median filter diff of %f pct' %
                  (var, kernel_size, ds.attrs[var + '_med_diff_pct']))
            filtered = scipy.signal.medfilt(ds[var], kernel_size=kernel_size)
            bads = (100 * np.abs(ds[var] - filtered)/ds[var] >
                    ds.attrs[var + '_med_diff_pct'])
            ds[var][bads] = np.nan

            notetxt = ('Values filled where percent difference between '
                       '%d-point median filter and original values is greater '
                       'than %f. ' %
                       (kernel_size, ds.attrs[var + '_med_diff_pct']))

            if 'note' in ds[var].attrs:
                ds[var].attrs['note'] = notetxt + ds[var].attrs['note']
            else:
                ds[var].attrs.update({'note': notetxt})

        if var + '_bad_ens' in ds.attrs:
            inc = np.arange(0, len(ds.attrs[var + '_bad_ens']), 2)
            for n in inc:
                print('%s: Trimming using bad_ens %d:%d' %
                      (var,
                       ds.attrs[var + '_bad_ens'][n],
                       ds.attrs[var + '_bad_ens'][n+1]))
                bads = np.full(ds[var].shape, False)
                bads[np.arange(ds.attrs[var + '_bad_ens'][n],
                               ds.attrs[var + '_bad_ens'][n+1])] = True
                ds[var][bads] = np.nan

    return ds
