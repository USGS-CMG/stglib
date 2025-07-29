Sea-Bird instruments
********************

SBE 37 MicroCAT
===============

Starting from exported .asc file, use :doc:`runots.py </runots>` to process using the two :doc:`configuration files </config>`.


SBE 26plus Seagauge (tides)
===========================

Starting from exported .tid or .wb file, use :doc:`runots.py </runots>` to process using the two :doc:`configuration files </config>`.

In the cdf2nc step, you can optionally input a .nc file that contains salinity and water temperature to be used to improve the estimates of water_level and water_depth variables by including [--salwtemp SALWTEMP]. Note: SALWTEMP nc file needs to include time period of the deployment being processed.



SBE 26plus Seagauge (waves)
===========================

Starting from exported .wb file, use :doc:`runots.py </runots>` to process using the two :doc:`configuration files </config>`.
