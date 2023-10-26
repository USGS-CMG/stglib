# stglib - Process data from a variety of oceanographic instrumentation

[![Documentation Status](https://readthedocs.org/projects/stglib/badge/?version=latest)](http://stglib.readthedocs.io/en/latest/?badge=latest)
![stglib](https://github.com/USGS-CMG/stglib/workflows/stglib/badge.svg)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/stglib/badges/version.svg)](https://anaconda.org/conda-forge/stglib)

This package contains code to process data from a variety of oceanographic instrumentation, consistent with the procedures of the USGS [Coastal/Marine Hazards and Resources Program](https://marine.usgs.gov) (formerly Coastal and Marine Geology Program).

Currently, this package has at least partial support for:

- Nortek Aquadopp profilers, in mean-current wave-burst, and HR modes
- Nortek Vector velocimeters
- Nortek Signature profilers
- RBR pressure (including waves) and turbidity sensors
- YSI EXO2 water-quality sondes
- SonTek IQ flow monitors
- WET labs sensors, including ECO NTUSB and ECO PAR
- Onset HOBO pressure sensors
- Vaisala Weather Transmitter WXT sensors
- In-Situ Aqua TROLL sensors
- RD Instruments ADCPs
- Moving-boat ADCP data processed using [QRev](https://hydroacoustics.usgs.gov/movingboat/QRev.shtml), for use in index-velocity computation
- EofE ECHOLOGGER altimeters

We have plans to support:

- RDI Sentinel V profilers

This package makes heavy use of [NumPy](http://www.numpy.org), [xarray](http://xarray.pydata.org/en/stable/), and [netCDF4](http://unidata.github.io/netcdf4-python/). It works on Python 3.9+.

[Read the documentation](http://stglib.readthedocs.io/).
