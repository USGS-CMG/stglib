Onset HOBO
**********

Starting from the exported .csv file from HOBOware software, use :doc:`runots.py </runots>` to process using the two :doc:`configuration files </config>`.

In the cdf2nc step, you can optionally input a .nc file that contains salinity and water temperature to be used to improve the estimates of water_level variable by including [--salwtemp SALWTEMP]. Note: SALWTEMP nc file needs to include time period of the deployment being processed.

