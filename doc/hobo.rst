Onset HOBO
**********

Data will generally be processed using a series of run scripts. The first script for each instrument type
depends on two :doc:`configuration files </config>`.

Instrument data to raw .cdf
===========================

Convert from exported CSV file to a raw netCDF file with .cdf extension using ``runhwlbcsv2cdf.py``.

runhwlbcsv2cdf.py
-----------------

.. argparse::
   :ref: stglib.core.cmd.hwlbcsv2cdf_parser
   :prog: runhwlbcsv2cdf.py

Raw .cdf to CF-compliant .nc
============================

Convert the raw .cdf data into CF-compliant netCDF file with .nc extension using ``runhwlbcdf2nc.py``.

runhwlbcdf2nc.py
----------------

.. argparse::
   :ref: stglib.core.cmd.hwlbcdf2nc_parser
   :prog: runhwlbcdf2nc.py
