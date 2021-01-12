import pandas as pd
import xarray as xr
import matplotlib.dates

from ..core import utils
from . import rdradcp


def raw_to_cdf(metadata):
    """Load a Aquadopp text files and output to netCDF format"""

    # TODO: clock drift code
    # TODO: logmeta code

    basefile = metadata["basefile"]

    if "prefix" in metadata:
        basefile = metadata["prefix"] + basefile

    # utils.check_valid_metadata(metadata)

    adcp = rdradcp.rdradcp(basefile, 1)

    # return adcp

    ds = xr.Dataset()

    # write out metadata first, then deal exclusively with xarray attrs
    ds = utils.write_metadata(ds, metadata)

    del metadata

    ds['time'] = xr.DataArray(pd.DatetimeIndex(matplotlib.dates.num2date(adcp.mtime)).tz_convert(None), dims='time')
    ds['cell'] = xr.DataArray(range(adcp.config.n_cells), dims='cell')
    ds['beam'] = xr.DataArray(range(adcp.config.n_beams), dims='beam')

    # not used:
    # 'latitude', 'longitude', 'name', 'ensemble_data', 'bin_data',
    # 'bt_ampl', 'bt_corr', 'bt_mode', 'bt_perc_good', 'bt_range', 'bt_vel'

    for var in ['corr', 'intens', 'perc_good', 'status']:
        ds[var] = xr.DataArray(getattr(adcp, var), dims=['time', 'cell', 'beam'])

    # rdradcp returns 'east_vel' etc even if it's in beam coordinates
    if adcp.config.coord_sys == 'beam':
        velvars = {'east_vel': 'vel_1', 'north_vel': 'vel_2', 'vert_vel': 'vel_3', 'error_vel': 'vel_4'}
    else:
        velvars = {'east_vel': 'east_vel', 'north_vel': 'north_vel', 'vert_vel': 'vert_vel', 'error_vel': 'error_vel'}

    for var in velvars:
        ds[velvars[var]] = xr.DataArray(getattr(adcp, var), dims=['time', 'cell'])

    for var in ['heading', 'heading_std', 'pitch', 'pitch_std', 'roll', 'roll_std', 'temperature', 'pressure', 'pressure_std', 'depth', 'salinity', 'number']:
        ds[var] = xr.DataArray(getattr(adcp, var), dims='time')

    for k in dir(adcp.config):
        if '__' not in k:
            ds.attrs[k] = getattr(adcp.config, k)

    if utils.is_cf(ds):
        pass
    else:
        print("about to create epic times")
        ds = utils.create_epic_times(ds)

    # configure file
    if "prefix" in ds.attrs:
        cdf_filename = ds.attrs["prefix"] + ds.attrs["filename"] + "-raw.cdf"
    else:
        cdf_filename = ds.attrs["filename"] + "-raw.cdf"

    ds.to_netcdf(cdf_filename, unlimited_dims=["time"])

    print("Finished writing data to %s" % cdf_filename)

    return ds
