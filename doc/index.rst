.. stglib documentation master file, created by
   sphinx-quickstart on Mon Mar 19 15:05:48 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

stglib Documentation
====================

stglib is a library of routines to process data from a variety of oceanographic instrumentation, consistent with the procedures of the Sediment Transport Group at the USGS Woods Hole Coastal and Marine Science Center.

stglib serves two distinct but related purposes:

1. A library to import data into xarray Datasets, which can be manipulated in an interactive Python environment (e.g., IPython, Jupyter notebooks).
2. A series of run scripts to process data into EPIC-compliant netCDF files for release to the public.


Currently, this package has at least partial support for:

- Nortek Aquadopp profilers, in mean-current and wave-burst modes
- RBR d|wave pressure sensors
- YSI EXO2 water-quality sondes
- SonTek IQ flow monitors
- WET labs sensors, including ECO NTUSB and ECO PAR
- Onset HOBO pressure sensors
- Moving-boat ADCP data processed using QRev_, for use in index-velocity computation

.. _QRev: https://hydroacoustics.usgs.gov/movingboat/QRev.shtml

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   install
   overview
   config
   atmos
   aqd
   wvs
   rsk
   exo
   iq
   wet
   hobo
   indexvel
   turnaround
   code
   contributing


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
