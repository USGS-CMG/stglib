RBR instruments
***************

Data will generally be processed using a series of run scripts. The first script for each instrument type
depends on two :doc:`configuration files </config>`.

Instrument data to raw .cdf
===========================

First, export data to text (.zip) format (File -> Export -> Text (\*.zip)) from the Ruskin software.

Then use :doc:`runots.py </runots>` to process using the two :doc:`configuration files </config>`.

Waves processing
================

Option 1: stglib internal code
------------------------------

This option is recommended as it is significantly faster than using pyDIWASP option. Generate the waves statistics and incorporate them into an CF-compliant netCDF file with .nc extension using :doc:`runots.py </runots>`.

Option 2: pyDIWASP
------------------

Generate the waves statistics using pyDIWASP and incorporate them into an CF-compliant netCDF file with .nc extension using :doc:`runots.py </runots>`. This method is slower than option 1.
