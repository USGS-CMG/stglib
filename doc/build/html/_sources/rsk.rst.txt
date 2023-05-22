RBR instruments
***************

Data will generally be processed using a series of run scripts. The first script for each instrument type
depends on two :doc:`configuration files </config>`.

Instrument data to raw .cdf
===========================

First, export data to text (.zip) format (File -> Export -> Text (\*.zip)) from the Ruskin software.

Then use ``runrskcsv2cdf.py`` to convert this exported data to a raw netCDF .cdf file.

runrskcsv2cdf.py
----------------

.. argparse::
   :ref: stglib.core.cmd.rskcsv2cdf_parser
   :prog: runrskcsv2cdf.py

Raw .cdf to CF-compliant .nc
============================

Convert the raw .cdf data into an CF-compliant netCDF file with .nc extension using ``runrskcdf2nc.py``, optionally including :doc:`atmospheric correction </atmos>` of the pressure data.

runrskcdf2nc.py
---------------

.. argparse::
   :ref: stglib.core.cmd.rskcdf2nc_parser
   :prog: runrskcdf2nc.py


Waves processing
================

Option 1: internal waves code
-----------------------------

This option is recommended as it does not require MATLAB. Generate the waves statistics and incorporate them into an CF-compliant netCDF file with .nc extension using ``runrsknc2waves.py``.

runrsknc2waves.py
~~~~~~~~~~~~~~~~~

.. argparse::
  :ref: stglib.core.cmd.rsknc2waves_parser
  :prog: runrsknc2waves.py


Option 2: DIWASP
----------------

Run DIWASP (within MATLAB) to produce wave statistics (see ``scripts/rundiwasp.m`` for an example run script). DIWASP must be run within MATLAB.

Incorporate the DIWASP statistics into an CF-compliant netCDF file with .nc extension using ``runrsknc2diwasp.py``.

runrsknc2diwasp.py
~~~~~~~~~~~~~~~~~~
