# stglib - Process data from a variety of oceanographic instrumentation

This package contains code to process data from a variety of oceanographic instrumentation, consistent with the procedures of the Sediment Transport Group at the USGS Woods Hole Coastal and Marine Science Center.

Currently, this package has at least partial support for:

- Nortek Aquadopp profilers, in mean-current and wave-burst modes
- RBR d|wave pressure sensors
- YSI EXO2 water-quality sondes
- SonTek IQ flow monitors
- WET labs sensors, including ECO NTUSB and ECO PAR
- Onset HOBO pressure sensors

This package makes heavy use of [NumPy](http://www.numpy.org), [xarray](http://xarray.pydata.org/en/stable/), and [netCDF4](http://unidata.github.io/netcdf4-python/). It has been tested only on Python 3.6+.

To use this package, `import stglib` and see below for relevant usages.

## Nortek Aquadopp

### Mean-current mode

Processing consists of two main steps:

1. Convert from text to a raw netCDF file with `.cdf` extension (`scripts/runaqdhdr2cdf.py`)
2. Convert the raw `.cdf` data into an EPIC-compliant netCDF file with `.nc` extension (`scripts/runaqdcdf2nc.py`), optionally including atmospheric correction of the pressure data (see `scripts/aqd_make_press_ac.ipynb`)

### Wave-burst mode

Processing consists of three main steps:

1. Convert from text to a raw netCDF file with `.cdf` extension (`scripts/runwvswad2cdf.py`)
2. Convert the raw `.cdf` data into an EPIC-compliant netCDF file with `.nc` extension (`scripts/runwvscdf2nc.py`), optionally including atmospheric correction of the pressure data (see `scripts/aqd_make_press_ac.ipynb`)
3. Run DIWASP (within MATLAB) to produce wave statistics (see `scripts/rundiwasp_aqd.m`), and incorporate these statistics into an EPIC-compliant netCDF file with `.nc` extension (`scripts/runwvsnc2diwasp.py`). Note that DIWASP is a MATLAB package and must be run from MATLAB before using this module. This code uses a version of DIWASP that has been updated to use radial beam velocities in the computation of wave statistics.

## RBR d|wave

Processing consists of three main steps:

1. Convert from `.rsk` binary to a raw netCDF file with `.cdf` extension (`scripts/runrskrsk2cdf.py`)
2. Convert the raw `.cdf` data into an EPIC-compliant netCDF file with `.nc` extension (`scripts/runrskcdf2nc.py`), optionally including atmospheric correction of the pressure data (see `scripts/dw_make_press_ac.ipynb`)
3. Run DIWASP (within MATLAB) to produce wave statistics (see `scripts/rundiwasp.m`), and incorporate these statistics into an EPIC-compliant netCDF file with `.nc` extension (`scripts/runrsknc2diwasp.py`). Note that DIWASP is a MATLAB package and must be run from MATLAB before using this module. A sample MATLAB run file for DIWASP is included in the `scripts` directory.

## YSI EXO2

Currently this module supports reading the `.csv` file exported from the KOR software into an xarray `Dataset` using `stglib.exo.read_exo()`

## SonTek IQ

Currently this module supports reading the `.mat` file exported from the SonTek IQ software into an xarray `Dataset` using `stglib.iq.read_iq()`

## WET labs ECO sensors

### NTUSB

Currently this module supports reading the text file saved from the terminal program used to interface with the instrument into an xarray `Dataset` using `stglib.eco.read_ntu()`

### PAR

Currently this module supports reading the text file saved from the terminal program used to interface with the instrument into an xarray `Dataset` using `stglib.eco.read_par()`

## Onset HOBO

Currently this module supports reading the `.csv` file exported from the HOBO software into an xarray `Dataset` using `stglib.hobo.read_hobo()`
