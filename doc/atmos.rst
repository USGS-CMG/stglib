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
      comment: Atmospherically corrected using NOAA station 9414523 Redwood City, CA
  unlimited dimensions:
  current shape = (8874,)
  filling on

A second optional but highly recommended attribute for the ``atmpres`` variable is ``comment``. This can be used to store information about the station used for atmospheric compensation. For example, ``atmpres:comment = "Atmospherically corrected using NOAA station 9414523 Redwood City, CA"``.

The time base of the atmospheric pressure file must be the same as that of the instrument pressure record. This file will be used by the run script to atmospherically compensate the pressure record.

Steps to generate an atmospheric pressure file
----------------------------------------------

The Jupyter notebooks ``aqd_make_press_ac.ipynb``, ``dw_make_press_ac.ipynb``, and ``exo_make_press_ac.ipynb`` may help you with this process; the outline is described below.

1. Obtain timeseries of atmospheric pressure. Useful sources include `NWIS <https://nwis.waterdata.usgs.gov/nwis>`_ (parameter code 00025), `MesoWest <http://mesowest.utah.edu>`_, `NOAA Tides & Currents <https://tidesandcurrents.noaa.gov>`_, and `NERR sites <https://cdmo.baruch.sc.edu>`_.

2. Get the atmospheric pressure onto the same time base, with the same number of samples, as the instrument pressure record. For burst data, the atmospheric pressure should have the same length as the number of bursts.

3. Generate the atmospheric-compensated pressure record by subtracting atmospheric and a suitable offset so that the instrument reads as close to zero as possible in air (before/after deployment). That is, ``P_1ac = P_1 - atmos - offset``. You will probably iterate on values for the offset until you end up with a ``P_1ac`` that has near-zero values in air.

4. The resulting ``atmpres.cdf`` file (or any name of your choosing) should use CF conventions for time and should match the structure as reported by ``ncinfo`` shown above. Routines within stglib will read the atmospheric data and perform the compensation as part of the processing.
