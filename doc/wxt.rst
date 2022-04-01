Processing Vaisala WXT536 data
******************************

Data will generally be processed using a series of run scripts. The first script for each instrument type
depends on two :doc:`configuration files </config>`.

Instrument data to raw .cdf
===========================

Convert from exported .csv file to a raw netCDF file with .cdf extension using ``runwxtcsv2cdf.py``.

runwxtcsv2cdf.py
----------------

.. argparse::
   :ref: stglib.core.cmd.wxtcsv2cdf_parser
   :prog: runwxtcsv2cdf.py

Raw .cdf to clean .nc
=====================

Convert the raw .cdf data into a CF-compliant netCDF file with .nc extension using ``runwxtcdf2nc.py``.

runwxtcdf2nc.py
---------------

.. argparse::
   :ref: stglib.core.cmd.wxtcdf2nc_parser
   :prog: runwxtcdf2nc.py
