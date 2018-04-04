Atmospheric compensation
========================

Pressure can be atmospherically compensated using a special .cdf file in each instrument's directory. Atmospheric pressure data is usually obtained from a nearby weather station, or a sensor deployed specifically for this purpose as part of the study. The .cdf file should have the following general structure:

::

  <class 'netCDF4._netCDF4.Dataset'>
  root group (NETCDF4 data model, file format HDF5):
    dimensions(sizes): time(8874)
    variables(dimensions): int64 time(time), float64 atmpres(time)
    groups:

The file contains one dimension, ``time``, and one variable, ``atmpres``. The ``atmpres`` variable has one attribute, ``offset``,
which sets the zero offset for the pressure sensor. In the example below, the offset is -10.25 dbar.

::

  <class 'netCDF4._netCDF4.Variable'>
  float64 atmpres(time)
      _FillValue: nan
      offset: -10.25
  unlimited dimensions:
  current shape = (8874,)
  filling on

The time base of the atmospheric pressure file must be the same as that of the instrument pressure record. This file will be used by the run scripts to atmospherically compensate the pressure record.
