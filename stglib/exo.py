from __future__ import division, print_function
import warnings
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

    try:
        exo = pd.read_csv(filnam,
                          skiprows=skiprows,
                          infer_datetime_format=True,
                          parse_dates=[['Date (MM/DD/YYYY)',
                                        'Time (HH:MM:SS)']],
                          encoding=encoding)
    except UnicodeDecodeError:
        exo = pd.read_csv(filnam,
                          skiprows=skiprows,
                          infer_datetime_format=True,
                          parse_dates=[['Date (MM/DD/YYYY)',
                                        'Time (HH:MM:SS)']],
                          encoding='mac-roman')
    except ValueError:
        print((' *** Could not decode header. '
               'Have you specified skiprows correctly?'))

    exo.rename(columns={'Date (MM/DD/YYYY)_Time (HH:MM:SS)': 'time'},
               inplace=True)
    exo.set_index('time', inplace=True)
    exo.rename(columns=lambda x: x.replace(' ', '_'), inplace=True)
    exo.rename(columns=lambda x: x.replace('/', '_per_'), inplace=True)
    exo['Press_dbar'] = exo['Press_psi_a'] * 0.689476

    exo = xr.Dataset(exo)
    hdr = read_exo_header(filnam, encoding=encoding)
    exo.attrs['serial_number'] = hdr['serial_number']
    exo.attrs['INST_TYPE'] = 'YSI EXO2 Multiparameter Sonde'
    exo.attrs['COMPOSITE'] = 0

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

    metadata.pop('skiprows')

    # write out metadata first, then deal exclusively with xarray attrs
    ds = utils.write_metadata(ds, metadata)

    del metadata

    ds = utils.create_epic_times(ds)

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
    ds = xr.open_dataset(cdf_filename).load()
    ds.close()

    # Clip data to in/out water times or via good_ens
    ds = utils.clip_ds(ds)

    ds = ds_rename_vars(ds)
    if 'C_51' in ds:
        ds['C_51'].values = ds['C_51'].values/10  # convert from mS/cm to S/m

    # ds = ds_add_attrs(ds)

    ds = ds.drop(['Press_psi_a',
                  'Site_Name',
                  'Fault_Code',
                  'Time_(Fract._Sec)'])

    if atmpres:
        print("Atmospherically correcting data")

        met = xr.open_dataset(atmpres).load()
        met.close()
        # need to save attrs before the subtraction, otherwise they are lost
        attrs = ds['P_1'].attrs
        ds['P_1ac'] = ds['P_1'] - met['atmpres'] - met['atmpres'].offset
        print('Correcting using offset of %f' % met['atmpres'].offset)
        ds['P_1ac'].attrs = attrs

    ds = exo_qaqc(ds)

    # assign min/max:
    ds = utils.add_min_max(ds)

    ds = utils.create_epic_times(ds)

    ds = exo_add_delta_t(ds)

    # add lat/lon coordinates
    ds = utils.ds_add_lat_lon(ds)

    ds = ds_add_attrs(ds)

    # ds = utils.create_water_depth(ds)
    ds = utils.create_nominal_instrument_depth(ds)

    ds = exo_create_depth(ds)

    # add lat/lon coordinates to each variable
    for var in ds.variables:
        if (var not in ds.coords) and ('time' not in var):
            ds = utils.add_lat_lon(ds, var)
            if 'P_1' not in var:
                ds = exo_add_depth(ds, var)
            # cast as float32
            ds = utils.set_var_dtype(ds, var)

    ds = utils.rename_time(ds)

    # No longer report depth
    ds = ds.drop('Depth_m')

    # Write to .nc file
    print("Writing cleaned/trimmed data to .nc file")
    nc_filename = ds.attrs['filename'] + '-a.nc'

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
                'Temp_∞C': 'T_28',
                'Cond_mS_per_cm': 'C_51',
                'SpCond_mS_per_cm': 'SpC_48',
                'Sal_psu': 'S_41',
                'ODO_%_sat': 'OST_62',
                'ODO_mg_per_L': 'DO',
                'Turbidity_NTU': 'Turb',
                'pH': 'pH_159',
                'pH_mV': 'pHmV'}

    # check to make sure they exist before trying to rename
    newvars = {}
    for k in varnames:
        if k in ds:
            newvars[k] = varnames[k]

    return ds.rename(newvars)


def ds_add_attrs(ds):
    # Update attributes for EPIC and STG compliance
    ds = utils.ds_coord_no_fillvalue(ds)

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
        if (var not in ds.coords) and ('time' not in var):
            add_attributes(ds[var], ds.attrs)

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

    # S_41 needs to be first in list for trim_by_salinity()
    for var in ['S_41', 'C_51', 'SpC_48', 'T_28', 'Turb', 'fDOMRFU', 'fDOMQSU', 'CHLrfu', 'Fch_906', 'BGAPErfu', 'BGAPE', 'OST_62', 'DO', 'pH_159', 'pHmV']:
        ds = trim_min(ds, var)

        ds = trim_max(ds, var)

        ds = trim_min_diff(ds, var)

        ds = trim_max_diff(ds, var)

        ds = trim_med_diff(ds, var)

        ds = trim_med_diff_pct(ds, var)

        ds = trim_bad_ens(ds, var)

        ds = trim_by_salinity(ds, var) # this must come last

    return ds


def trim_min(ds, var):
    if var + '_min' in ds.attrs:
        print('%s: Trimming using minimum value of %f' %
              (var, ds.attrs[var + '_min']))
        ds[var][ds[var] < ds.attrs[var + '_min']] = np.nan

        notetxt = ('Values filled where less than %f units. ' %
                   ds.attrs[var + '_min'])

        ds = insert_note(ds, var, notetxt)

    return ds


def trim_max(ds, var):
    if var + '_max' in ds.attrs:
        print('%s: Trimming using maximum value of %f' %
              (var, ds.attrs[var + '_max']))
        ds[var][ds[var] > ds.attrs[var + '_max']] = np.nan

        notetxt = ('Values filled where greater than %f units. ' %
                   ds.attrs[var + '_max'])

        ds = insert_note(ds, var, notetxt)

    return ds


def trim_min_diff(ds, var):
    if var + '_min_diff' in ds.attrs:
        print('%s: Trimming using minimum diff of %f' %
              (var, ds.attrs[var + '_min_diff']))
        ds[var][np.ediff1d(
            ds[var], to_begin=0) < ds.attrs[var + '_min_diff']] = np.nan

        notetxt = ('Values filled where data decreases by more than %f '
                   'units in a single time step. ' %
                   ds.attrs[var + '_min_diff'])

        ds = insert_note(ds, var, notetxt)

    return ds


def trim_max_diff(ds, var):
    if var + '_max_diff' in ds.attrs:
        print('%s: Trimming using maximum diff of %f' %
              (var, ds.attrs[var + '_max_diff']))
        ds[var][np.ediff1d(
            ds[var], to_begin=0) > ds.attrs[var + '_max_diff']] = np.nan

        notetxt = ('Values filled where data increases by more than %f '
                   'units in a single time step. ' %
                   ds.attrs[var + '_max_diff'])

        ds = insert_note(ds, var, notetxt)

    return ds


def trim_med_diff(ds, var):
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

        ds = insert_note(ds, var, notetxt)

    return ds


def trim_med_diff_pct(ds, var):
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

        ds = insert_note(ds, var, notetxt)

    return ds


def trim_bad_ens(ds, var):
    if var + '_bad_ens' in ds.attrs:
        inc = np.arange(0, len(ds.attrs[var + '_bad_ens']), 2)
        for n in inc:
            print('%s: Trimming using bad_ens %s' %
                  (var, str(ds.attrs[var + '_bad_ens'][n:n+2])))
            if isinstance(ds.attrs[var + '_bad_ens'][n], str):
                bads = ((ds['time'] >= np.datetime64(ds.attrs[var + '_bad_ens'][n])) &
                        (ds['time'] <= np.datetime64(ds.attrs[var + '_bad_ens'][n+1])))
                ds[var] = ds[var].where(~bads)
            else:
                bads = np.full(ds[var].shape, False)
                bads[np.arange(ds.attrs[var + '_bad_ens'][n],
                               ds.attrs[var + '_bad_ens'][n+1])] = True
                ds[var][bads] = np.nan

            notetxt = "Data clipped using bad_ens values of %s. " % (
                str(ds.attrs[var + '_bad_ens']))

            ds = insert_note(ds, var, notetxt)

    return ds


def trim_by_salinity(ds, var):
    if 'trim_by_salinity' in ds.attrs and ds.attrs['trim_by_salinity'].lower() == 'true': # xarray doesn't support writing attributes as booleans
        print('%s: Trimming using valid salinity threshold' % var)
        ds[var][ds['S_41'].isnull()] = np.nan

        notetxt = 'Values filled using valid salinity threshold. '

        ds = insert_note(ds, var, notetxt)

    return ds


def insert_note(ds, var, notetxt):
    if 'note' in ds[var].attrs:
        ds[var].attrs['note'] = notetxt + ds[var].attrs['note']
    else:
        ds[var].attrs.update({'note': notetxt})

    return ds


def exo_create_depth(ds):

    depth = ds.attrs['WATER_DEPTH'] - ds.attrs['initial_instrument_height']
    ds['depth'] = xr.DataArray([depth], dims='depth')
    ds['depth'].attrs['positive'] = 'down'
    ds['depth'].attrs['axis'] = 'z'
    ds['depth'].attrs['units'] = 'm'
    ds['depth'].attrs['epic_code'] = 3
    ds['depth'].encoding['_FillValue'] = 1e35

    return ds


def exo_add_depth(ds, var):
    ds[var] = xr.concat([ds[var]], dim=ds['depth'])

    # Reorder so lat, lon are at the end.
    dims = [d for d in ds[var].dims if (d != 'depth')]
    dims.extend(['depth'])
    dims = tuple(dims)

    ds[var] = ds[var].transpose(*dims)

    return ds


def exo_add_delta_t(ds):
    deltat = np.asscalar(
        (ds['time'][1] - ds['time'][0]) / np.timedelta64(1, 's'))
    if not deltat.is_integer():
        warnings.warn('DELTA_T is not an integer; casting as int in attrs')

    ds.attrs['DELTA_T'] = int(deltat)

    return ds
