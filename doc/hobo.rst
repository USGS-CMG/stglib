Onset HOBO
**********

Data will generally be processed using a series of run scripts. The first script for each instrument type
depends on two :doc:`configuration files </config>`.

Instrument data to raw .cdf
===========================

Convert from exported CSV file to a raw netCDF file with .cdf extension using ``runhobocsv2cdf.py``.

runhobocsv2cdf.py
-----------------

.. argparse::
   :ref: stglib.core.cmd.hobocsv2cdf_parser
   :prog: runhobocsv2cdf.py

Raw .cdf to CF-compliant .nc
============================

Convert the raw .cdf data into CF-compliant netCDF file with .nc extension using ``runhobocdf2nc.py``.

runhobocdf2nc.py
----------------

.. argparse::
   :ref: stglib.core.cmd.hobocdf2nc_parser
   :prog: runhobocdf2nc.py
