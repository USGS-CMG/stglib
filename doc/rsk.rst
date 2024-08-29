RBR instruments
***************

Data will generally be processed using a series of run scripts. The first script for each instrument type
depends on two :doc:`configuration files </config>`.

Instrument data to raw .cdf
===========================

First, export data to text (.zip) format (File -> Export -> Text (\*.zip)) from the Ruskin software.

Then use :doc:`runots </runots>` to process using the two :doc:`configuration files </config>`.

Waves processing
================

Option 1: internal waves code
-----------------------------

This option is recommended as it does not require MATLAB. Generate the waves statistics and incorporate them into an CF-compliant netCDF file with .nc extension using :doc:`runots </runots>`.

Option 2: DIWASP
----------------

Run DIWASP (within MATLAB) to produce wave statistics (see ``scripts/rundiwasp.m`` for an example run script). DIWASP must be run within MATLAB.

Incorporate the DIWASP statistics into an CF-compliant netCDF file with .nc extension using ``runrsknc2diwasp.py``.

runrsknc2diwasp.py
~~~~~~~~~~~~~~~~~~
