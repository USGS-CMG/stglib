# stglib - Process data from a variety of oceanographic instrumentation

This package contains code to process data from a variety of oceanographic instrumentation, consistent with the procedures of the Sediment Transport Group at the USGS Woods Hole Coastal and Marine Science Center.

Currently, this package has at least partial support for:

- Nortek Aquadopp profilers, in mean-current and wave-burst modes
- RBR d|wave pressure sensors
- YSI EXO2 water-quality sondes
- SonTek IQ flow monitors
- WET labs ECO NTUSB turbidity sensors
- Onset HOBO pressure sensors

This package makes heavy use of [NumPy](http://www.numpy.org), [xarray](http://xarray.pydata.org/en/stable/), and [netCDF4](http://unidata.github.io/netcdf4-python/).

## Nortek Aquadopp

Processing consists of two main steps:

1. Convert from text to a raw netCDF file with `.cdf` extension
2. Convert the raw `.cdf` data into an EPIC-compliant netCDF file with `.nc` extension, optionally including atmospheric correction of the pressure data

### Text to raw netCDF (.cdf)

This step will generally be completed by using the import statement `from aqdlib import aqdhdr2cdf` and calling `aqdhdr2cdf.hdr_to_cdf()`, or by running `aqdhdr2cdf.py` from the command line.

### Raw netCDF (.cdf) to EPIC-compliant and processed netCDF (.nc)

This step will generally be completed by using the import statement `from aqdlib import aqdcdf2nc` and calling `aqdcdf2nc.cdf_to_nc()`, or by running `aqdcdf2nc.py` from the command line. When calling `cdf_to_nc()`, the user may provide the path to a netCDF file consisting of atmospheric pressure, which will be used to atmospherically correct the pressure data. This path can also be passed as a command-line argument to `aqdcdf2nc.py`.

## RBR d|wave

Processing consists of three main steps:

1. Convert from `.rsk` binary to a raw netCDF file with `.cdf` extension
2. Convert the raw `.cdf` data into an EPIC-compliant netCDF file with `.nc` extension, optionally including atmospheric correction of the pressure data
3. Run DIWASP (within MATLAB) to produce wave statistics, and incorporate these statistics into an EPIC-compliant netCDF file with `.nc` extension

### Raw binary to raw netCDF (.cdf)

This step will generally be completed by using the import statement `from rsklib import rskrsk2cdf` and calling `rskrsk2cdf.hdr_to_cdf()`, or by running `rskrsk2cdf.py` from the command line.

### Raw netCDF (.cdf) to EPIC-compliant and processed netCDF (.nc)

This step will generally be completed by using the import statement `from rsklib import rskcdf2nc` and calling `rskcdf2nc.cdf_to_nc()`, or by running `rskcdf2nc.py` from the command line. When calling `cdf_to_nc()`, the user may provide the path to a netCDF file consisting of atmospheric pressure, which will be used to atmospherically correct the pressure data. This path can also be passed as a command-line argument to `rskcdf2nc.py`.

### DIWASP processing and creation of EPIC-compliant wave statistics netCDF (.nc)

This step will generally be completed by using the import statement `from rsklib import rsknc2diwasp` and calling `rsknc2diwasp.nc_to_diwasp()`, or by running `rsknc2diwasp.py` from the command line. Note that DIWASP is a MATLAB package and must be run from MATLAB before using this module. A sample MATLAB run file for DIWASP is included in the `scripts` directory.

## YSI EXO2

Currently this module supports reading the `.xlsx` file exported from the KOR software into an xarray Dataset.

## SonTek IQ

Currently this module supports reading the `.mat` file exported from the SonTek IQ software into an xarray Dataset.

## WET labs ECO NTUSB

Currently this module supports reading the text file saved from the terminal program used to interface with the instrument into an xarray Dataset.

## Onset HOBO

Currently this module supports reading the `.csv` file exported from the HOBO software into an xarray Dataset.
