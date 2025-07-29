RBR instruments
***************

Data will generally be processed using a series of run scripts. The first script for each instrument type
depends on two :doc:`configuration files </config>`.


Instrument data to raw .cdf
===========================

First, export data to text (.zip) format (File -> Export -> Text (\*.zip)) from the Ruskin software.

Then use :doc:`runots.py </runots>` to process using the two :doc:`configuration files </config>`.

In the cdf2nc step, you can optionally input a .nc file that contains salinity and water temperature to be used to improve the estimates of water_level and water_depth variables by including [--salwtemp SALWTEMP]. Note: SALWTEMP nc file needs to include time period of the deployment being processed.


Waves processing
================

Option 1: stglib internal code
------------------------------

This option is recommended as it is significantly faster than using pyDIWASP option. Generate the waves statistics and incorporate them into an CF-compliant netCDF file with .nc extension using :doc:`runots.py </runots>`.

Option 2: pyDIWASP
------------------

Generate the waves statistics using pyDIWASP and incorporate them into an CF-compliant netCDF file with .nc extension using :doc:`runots.py </runots>`. This method is slower than option 1.
