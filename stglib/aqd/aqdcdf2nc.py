#!/usr/bin/env python

from __future__ import division, print_function
import xarray as xr
import sys
sys.path.insert(0, '/Users/dnowacki/Documents/aqdlib')
import aqdlib.qaqc as qaqc
import netCDF4


def cdf_to_nc(cdf_filename, metadata, atmpres=False):
    """
    Load a "raw" .cdf file and generate a processed .nc file
    """

    # Load raw .cdf data
    VEL = load_cdf(cdf_filename, metadata, atmpres=atmpres)

    # Clip data to in/out water times or via good_ens
    VEL = clip_ds(VEL, metadata)

    # Create water_depth variables
    VEL, metadata = qaqc.create_water_depth(VEL, metadata)

    # Create depth variable depending on orientation
    VEL, T = qaqc.set_orientation(VEL, VEL['TransMatrix'].values, metadata)

    # Transform coordinates from, most likely, BEAM to ENU
    u, v, w = qaqc.coord_transform(VEL['VEL1'].values, VEL['VEL2'].values, VEL['VEL3'].values,
        VEL['Heading'].values, VEL['Pitch'].values, VEL['Roll'].values, T, VEL.attrs['AQDCoordinateSystem'])

    VEL['U'] = xr.DataArray(u, dims=('time', 'bindist'))
    VEL['V'] = xr.DataArray(v, dims=('time', 'bindist'))
    VEL['W'] = xr.DataArray(w, dims=('time', 'bindist'))

    VEL = qaqc.magvar_correct(VEL, metadata)

    VEL['AGC'] = (VEL['AMP1'] + VEL['AMP2'] + VEL['AMP3']) / 3

    VEL = qaqc.trim_vel(VEL, metadata)

    VEL = qaqc.make_bin_depth(VEL, metadata)

    # Reshape and associate dimensions with lat/lon
    for var in ['U', 'V', 'W', 'AGC', 'Pressure', 'Temperature', 'Heading', 'Pitch', 'Roll']:
        VEL = da_reshape(VEL, var)

    # swap_dims from bindist to depth
    VEL = ds_swap_dims(VEL)

    VEL = ds_rename(VEL)

    VEL = ds_drop(VEL)

    VEL = ds_add_attrs(VEL, metadata)

    # TODO: Need to add all global attributes from CDF to NC file (or similar)
    VEL = qaqc.add_min_max(VEL)

    VEL = qaqc.add_final_metadata(VEL)

    nc_filename = metadata['filename'] + '.nc'

    VEL.to_netcdf(nc_filename, unlimited_dims='time')
    print('Done writing netCDF file', nc_filename)

    # rename time variables after the fact to conform with EPIC/CMG standards
    rename_time(nc_filename)

    print('Renamed dimensions')

    return VEL


def load_cdf(cdf_filename, metadata, atmpres=False):
    """
    Load raw .cdf file and, optionally, an atmospheric pressure .cdf file
    """

    ds = xr.open_dataset(cdf_filename, autoclose=True)

    if atmpres is not False:
        p = xr.open_dataset(atmpres, autoclose=True)
        # TODO: check to make sure this data looks OK
        # need to call load for waves; it's not in memory and throws error
        ds['Pressure_ac'] = xr.DataArray(ds['Pressure'].load() - (p['atmpres'] - p['atmpres'].offset))

    return ds


def clip_ds(ds, metadata):
    """
    Clip an xarray Dataset from metadata, either via good_ens or
    Deployment_date and Recovery_date
    """

    print('first burst in full file:', ds['time'].min().values)
    print('last burst in full file:', ds['time'].max().values)

    # clip either by ensemble indices or by the deployment and recovery date specified in metadata
    if 'good_ens' in metadata:
        # we have good ensemble indices in the metadata
        print('Clipping data using good_ens')

        ds = ds.isel(time=slice(metadata['good_ens'][0], metadata['good_ens'][1]))

        histtext = 'Data clipped using good_ens values of ' + metadata['good_ens'][0] + ', ' + metadata['good_ens'][1] + '. '
        if 'history' in ds.attrs:
            ds.attrs['history'] = histtext + ds.attrs['history']
        else:
            ds.attrs['history'] = histtext

    elif 'Deployment_date' in metadata and 'Recovery_date' in metadata:
        # we clip by the times in/out of water as specified in the metadata
        print('Clipping data using Deployment_date and Recovery_date')

        ds = ds.sel(time=slice(metadata['Deployment_date'], metadata['Recovery_date']))

        histtext = 'Data clipped using Deployment_date and Recovery_date of ' + metadata['Deployment_date'] + ', ' + metadata['Recovery_date'] + '. '
        if 'history' in ds.attrs:
            ds.attrs['history'] = histtext + ds.attrs['history']
        else:
            ds.attrs['history'] = histtext
    else:
        # do nothing
        print('Did not clip data; no values specified in metadata')

    print('first burst in trimmed file:', ds['time'].min().values)
    print('last burst in trimmed file:', ds['time'].max().values)

    return ds


# TODO: add analog input variables (OBS, NTU, etc)


def rename_time(nc_filename):
    """
    Rename time variables. Need to use netCDF4 module since xarray seems to have
    issues with the naming of time variables/dimensions
    """

    nc = netCDF4.Dataset(nc_filename, 'r+')
    timebak = nc['epic_time'][:]
    nc.renameVariable('time', 'time_cf')
    nc.renameVariable('epic_time', 'time')
    nc.renameVariable('epic_time2', 'time2')
    nc.close()

    # need to do this in two steps after renaming the variable
    # not sure why, but it works this way
    nc = netCDF4.Dataset(nc_filename, 'r+')
    nc['time'][:] = timebak
    nc.close()


def ds_swap_dims(ds):

    ds = ds.swap_dims({'bindist': 'depth'})

    # need to swap dims and then reassign bindist to be a normal variable (no longer a coordinate)
    valbak = ds['bindist'].values
    ds = ds.drop('bindist')
    ds['bindist'] = xr.DataArray(valbak, dims='depth')

    return ds


def da_reshape(ds, var, waves=False):
    """
    Add lon and lat dimensions to DataArrays and reshape to conform to our
    standard order
    """

    # Add the dimensions using concat
    ds[var] = xr.concat([ds[var]], dim=ds['lon'])
    ds[var] = xr.concat([ds[var]], dim=ds['lat'])

    # Reshape using transpose depending on shape
    if waves == False:
        if len(ds[var].shape) == 4:
            ds[var] = ds[var].transpose('time', 'lon', 'lat', 'bindist')
        elif len(ds[var].shape) == 3:
            ds[var] = ds[var].transpose('time', 'lon', 'lat')

    return ds


def ds_rename(ds, waves=False):
    """
    Rename DataArrays within Dataset for EPIC compliance
    """

    varnames = {'Pressure': 'P_1',
        'Temperature': 'Tx_1211',
        'Heading': 'Hdg_1215',
        'Pitch': 'Ptch_1216',
        'Roll': 'Roll_1217'}

    if 'Pressure_ac' in ds:
        varnames['Pressure_ac'] = 'P_1ac'

    if waves == False:
        varnames.update({'U': 'u_1205',
            'V': 'v_1206',
            'W': 'w_1204',
            'AGC': 'AGC_1202'})
    elif waves == True:
        varnames.update({'VEL1': 'vel1_1277',
            'VEL2': 'vel2_1278',
            'VEL3': 'vel3_1279',
            'AMP1': 'AGC1_1221',
            'AMP2': 'AGC2_1222',
            'AMP3': 'AGC3_1223'})

    ds.rename(varnames, inplace=True)

    return ds


def ds_drop(ds):
    """
    Drop old DataArrays from Dataset that won't make it into the final .nc file
    """

    todrop = ['VEL1',
        'VEL2',
        'VEL3',
        'AMP1',
        'AMP2',
        'AMP3',
        'Battery',
        'TransMatrix',
        'AnalogInput1',
        'AnalogInput2',
        'jd',
        'Depth']

    ds = ds.drop(todrop)

    return ds


def ds_add_attrs(ds, metadata, waves=False):
    """
    add EPIC and CMG attributes to xarray Dataset
    """

    def add_vel_attributes(vel, metadata):
        vel.attrs.update({'units': 'cm/s',
            'data_cmnt': 'Velocity in shallowest bin is often suspect and should be used with caution'})

        # TODO: why do we only do trim_method for Water Level SL?
        if 'trim_method' in metadata and metadata['trim_method'].lower() == 'water level sl':
            vel.attrs.update({'note': 'Velocity bins trimmed if out of water or if side lobes intersect sea surface'})

    def add_attributes(var, metadata, INFO):
        var.attrs.update({'serial_number': INFO['AQDSerial_Number'],
            'initial_instrument_height': metadata['initial_instrument_height'],
            'nominal_instrument_depth': metadata['nominal_instrument_depth'],
            'height_depth_units': 'm', 'sensor_type': INFO['INST_TYPE']})
        var.encoding['_FillValue'] = 1e35

    ds.attrs.update({'COMPOSITE': 0})

    # Update attributes for EPIC and STG compliance
    ds.lat.encoding['_FillValue'] = False
    ds.lon.encoding['_FillValue'] = False
    ds.depth.encoding['_FillValue'] = False
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

    ds['depth'].attrs.update({'units': 'm',
        'long_name': 'mean water depth',
        'initial_instrument_height': metadata['initial_instrument_height'],
        'nominal_instrument_depth': metadata['nominal_instrument_depth'],
        'epic_code': 3})

    if waves == False:
        ds['u_1205'].attrs.update({'name': 'u',
            'long_name': 'Eastward Velocity',
            'generic_name': 'u',
            'epic_code': 1205})

        ds['v_1206'].attrs.update({'name': 'v',
            'long_name': 'Northward Velocity',
            'generic_name': 'v',
            'epic_code': 1206})

        ds['w_1204'].attrs.update({'name': 'w',
            'long_name': 'Vertical Velocity',
            'generic_name': 'w',
            'epic_code': 1204})

        ds['AGC_1202'].attrs.update({'units': 'counts',
            'name': 'AGC',
            'long_name': 'Average Echo Intensity',
            'generic_name': 'AGC',
            'epic_code': 1202})

    elif waves == True:
        ds['vel1_1277'].attrs.update({'units': 'mm/s',
            'long_name': 'Beam 1 Velocity',
            'generic_name': 'vel1',
            'epic_code': 1277})

        ds['vel2_1278'].attrs.update({'units': 'mm/s',
            'long_name': 'Beam 2 Velocity',
            'generic_name': 'vel2',
            'epic_code': 1278})

        ds['vel3_1279'].attrs.update({'units': 'mm/s',
            'long_name': 'Beam 3 Velocity',
            'generic_name': 'vel3',
            'epic_code': 1279})

        ds['AGC1_1221'].attrs.update({'units': 'counts',
            'long_name': 'Echo Intensity (AGC) Beam 1',
            'generic_name': 'AGC1',
            'epic_code': 1221})

        ds['AGC2_1222'].attrs.update({'units': 'counts',
            'long_name': 'Echo Intensity (AGC) Beam 2',
            'generic_name': 'AGC2',
            'epic_code': 1222})

        ds['AGC3_1223'].attrs.update({'units': 'counts',
            'long_name': 'Echo Intensity (AGC) Beam 3',
            'generic_name': 'AGC3',
            'epic_code': 1223})

    ds['P_1'].attrs.update({'units': 'dbar',
        'name': 'P',
        'long_name': 'Pressure',
        'generic_name': 'depth',
        'epic_code': 1}) # TODO: is this generic name correct?

    if 'P_1ac' in ds:
        ds['P_1ac'].attrs.update({'units': 'dbar',
            'name': 'Pac',
            'long_name': 'Corrected pressure'})
        if 'P_1ac_note' in metadata:
            ds['P_1ac'].attrs.update({'note': metadata['P_1ac_note']})

        add_attributes(ds['P_1ac'], metadata, ds.attrs)

        ds.attrs['history'] = 'Atmospheric pressure compensated. ' + ds.attrs['history']

    ds['bin_depth'].attrs.update({'units': 'm',
        'name': 'bin depth'})

    if 'P_1ac' in ds:
        ds['bin_depth'].attrs.update({'note': 'Actual depth time series of velocity bins. Calculated as corrected pressure(P_1ac) - bindist.'})
    else:
        ds['bin_depth'].attrs.update({'note': 'Actual depth time series of velocity bins. Calculated as pressure(P_1) - bindist.'})

    ds['Tx_1211'].attrs.update({'units': 'C',
        'name': 'Tx',
        'long_name': 'Instrument Transducer Temperature',
        'generic_name': 'temp',
        'epic_code': 1211})

    ds['Hdg_1215'].attrs.update({'units': 'degrees',
        'name': 'Hdg',
        'long_name': 'Instrument Heading',
        'generic_name': 'hdg',
        'epic_code': 1215})

    if 'magnetic_variation_at_site' in metadata:
        ds['Hdg_1215'].attrs.update({'note': 'Heading is degrees true. Converted from magnetic with magnetic variation of ' + str(metadata['magnetic_variation_at_site'])})
    elif 'magnetic_variation' in metadata:
        ds['Hdg_1215'].attrs.update({'note': 'Heading is degrees true. Converted from magnetic with magnetic variation of ' + str(metadata['magnetic_variation'])})

    ds['Ptch_1216'].attrs.update({'units': 'degrees',
        'name': 'Ptch',
        'long_name': 'Instrument Pitch',
        'generic_name': 'ptch',
        'epic_code': 1216})

    ds['Roll_1217'].attrs.update({'units': 'degrees',
        'name': 'Roll',
        'long_name': 'Instrument Roll',
        'generic_name': 'roll',
        'epic_code': 1217})

    ds['bindist'].attrs.update({'units': 'm',
        'long_name': 'distance from transducer head',
        'blanking_distance': ds.attrs['AQDBlankingDistance'],
        'note': 'distance is along profile from instrument head to center of bin'})

    if waves == False:
        for v in ['AGC_1202', 'u_1205', 'v_1206', 'w_1204']:
            add_attributes(ds[v], metadata, ds.attrs)
        for v in ['u_1205', 'v_1206', 'w_1204']:
            add_vel_attributes(ds[v], metadata)
    elif waves == True:
        for v in ['vel1_1277', 'vel2_1278', 'vel3_1279', 'AGC1_1221', 'AGC2_1222', 'AGC3_1223']:
            add_attributes(ds[v], metadata, ds.attrs)

    for v in ['P_1', 'Tx_1211', 'Hdg_1215', 'Ptch_1216', 'Roll_1217', 'bin_depth', 'bindist']:
        add_attributes(ds[v], metadata, ds.attrs)

    return ds


def main():
    import sys
    sys.path.insert(0, '/Users/dnowacki/Documents/aqdlib')
    import aqdlib
    import argparse
    import yaml

    parser = argparse.ArgumentParser(description='Convert raw Aquadopp .cdf format to processed .nc files')
    parser.add_argument('cdfname', help='raw .CDF filename')
    parser.add_argument('gatts', help='path to global attributes file (gatts formatted)')
    parser.add_argument('config', help='path to ancillary config file (YAML formatted)')
    parser.add_argument('--atmpres', help='path to cdf file containing atmopsheric pressure data')

    args = parser.parse_args()

    # initialize metadata from the globalatts file
    metadata = aqdlib.read_globalatts(args.gatts)

    # Add additional metadata from metadata config file
    config = yaml.safe_load(open(args.config))

    for k in config:
        metadata[k] = config[k]

    if args.atmpres:
        VEL = cdf_to_nc(args.cdfname, metadata, atmpres=args.atmpres)
    else:
        VEL = cdf_to_nc(args.cdfname, metadata)

if __name__ == '__main__':
    main()
