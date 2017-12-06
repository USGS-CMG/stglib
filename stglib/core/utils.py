from __future__ import division, print_function
import csv

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

def write_metadata(ds, metadata):
    """Write out all metadata to CDF file"""

    for k in metadata:
        if k != 'instmeta': # don't want to write out instmeta dict, call it separately
            ds.attrs.update({k: metadata[k]})

    f = os.path.basename(inspect.stack()[1][1])

    ds.attrs.update({'history': 'Processed using ' + f + ' with Python ' +
        platform.python_version() + ', xarray ' + xr.__version__ + ', NumPy ' +
        np.__version__ + ', netCDF4 ' + netCDF4.__version__})

    return ds

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
