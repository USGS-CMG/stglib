# process the Virtuoso Tu from Matanzas 1110
import sys
sys.path.append('C:\\data\\Matanzas\\RBRVirtuosoTu')
#sys.path.append('/Users/dnowacki/Documents/python')
sys.path.append('C:\\projects\\python\\stglib')
sys.path.append('C:\\projects\\python\\xmltodict')
import stglib
from stglib.rsk.virtuosorsk2nc import virtuosorsk_to_cdf as rsk_to_cdf
from stglib.rsk.virtuosorsk2nc import virtuosocdf_to_nc as cdf_to_nc
#sys.path.append('C:\\projects\\python\\stglib\\rsk')
#from . import virtuosorsk2nc
import matplotlib.pyplot as plt
#from stglib import utils
import xarray as xr

gattfile = "..\\glob_att1110.txt"
metadata = stglib.read_globalatts(gattfile)
metadata['basefile'] = '054276_20180618_1614'
metadata['initial_instrument_height'] = 0.4 # best estimate from images
metadata['filename'] = '111012VTu' # this may change with new numbering

# make the raw data file
if 0:
    # ugly, until we get various dependencies settled
    xds, cdfname = rsk_to_cdf(metadata)
    #print(xds)
    
# amke the netcdf data file
if 0:
    cdf_to_nc(cdfname, writefile=True, format='netCDF4')

    # format='NETCDF3_64BIT' caused a problem
    # ValueError: could not safely cast array from dtype int64 to int32
    # in xarray netcdf3.py value = coerce_nc3_dtype(np.atleast_1d(value))
    
# do qaqc
if 1:
    # check for bad time base
    ds = xr.open_dataset(metadata['filename']+'-cal.nc', decode_times=False)
    buf = ds['time'].units.split()
    dt = ds['time'].diff(dim='time')
    str1 = 'Mean of the difference in time: \n%f %s' % (dt.mean(),buf[0])
    str2 = 'STD of the difference in time: \n%f %s' % (dt.std(),buf[0])
    
    fig, ax = plt.subplots(nrows=2,ncols=2,figsize=(10,10))
    ax[0,0].hist(dt,label='time difference')
    ax[0,0].legend()
    ax[0,0].set_xlabel(buf[0])
    ax[0,1].plot(dt)
    ax[0,1].set_xlabel('count')
    ax[0,1].set_ylabel(buf[0])
    ax[1,0].plot(ds['Turb'][:,0,0,0])
    ax[1,1].text(0.85, 0.5, str1,
        verticalalignment='bottom', horizontalalignment='right',
        transform=ax[1,1].transAxes,
        color='green', fontsize=10)
    ax[1,1].text(0.85, 0.3, str2,
        verticalalignment='bottom', horizontalalignment='right',
        transform=ax[1,1].transAxes,
        color='green', fontsize=10)
    fig.show()
        
    ds.close()

#fig, axes = plt.subplots(figsize=(15,2))
#line1 = axes.plot(dic['Turb'][:])
    
