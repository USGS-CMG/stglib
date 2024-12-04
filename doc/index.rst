.. stglib documentation master file, created by
   sphinx-quickstart on Mon Mar 19 15:05:48 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

stglib Documentation
====================

This package contains code to process data from a variety of oceanographic instrumentation, consistent with the procedures of the USGS Coastal/Marine Hazards and Resources Program (formerly Coastal and Marine Geology Program).

stglib serves two distinct but related purposes:

1. A library to import data into xarray Datasets, which can be manipulated in an interactive Python environment (e.g., IPython, Jupyter notebooks).
2. A series of run scripts to process data into CF-compliant netCDF files for release to the public.


Currently, this package has at least partial support for:

- Nortek Aquadopp profilers, in mean-current wave-burst, and HR modes
- Nortek Vector velocimeters
- Nortek Signature profilers
- Multiple RBR instruments, including moored and casting deployments
- YSI EXO2 water-quality sondes
- SonTek IQ flow monitors
- WET labs sensors, including ECO NTUSB and ECO PAR
- Onset HOBO pressure sensors
- Vaisala Weather Transmitter WXT sensors
- In-Situ Aqua TROLL sensors
- RD Instruments ADCPs
- Moving-boat ADCP data processed using QRev_, for use in index-velocity computation
- EofE ECHOLOGGER altimeters
- SBE 37 Microcat
- SBE 26plus Seagauge
- TruBlue pressure sensors

.. _QRev: https://hydroacoustics.usgs.gov/movingboat/QRev.shtml

We have plans to support:

- RDI Sentinel V profilers

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   install
   overview
   config
   Processing data with runots.py (use this!) <runots>
   atmos
   waves
   aqd
   aqdhr
   wvs
   rsk
   sig
   vec
   wet
   eofe
   exo
   hobo
   indexvel
   iq
   lisst
   mc
   sg
   tb
   tcm
   wxt
   turnaround
   code
   contributing


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
