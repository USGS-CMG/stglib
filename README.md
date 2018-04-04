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

[Read the documentation](http://stglib.readthedocs.io/).
