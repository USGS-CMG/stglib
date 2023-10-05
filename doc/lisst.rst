Sequoia Scientific LISST
************************

Data will generally be processed using a series of run scripts. The first script for each instrument type
depends on two :doc:`configuration files </config>`.

Exported ASCII data to raw .cdf
===============================

Convert from comma-delimited file to a raw netCDF file with .cdf extension using ``runlisstcsv2cdf.py``.

runlisstcsv2cdf.py
------------------

.. argparse::
   :ref: stglib.core.cmd.lisstcsv2cdf_parser
   :prog: runlisstcsv2cdf.py

Raw .cdf to CF-compliant .nc
============================

Convert the raw .cdf data into CF-compliant netCDF file with .nc extension using ``runlisstcdf2nc.py``.

runlisstcdf2nc.py
-----------------

.. argparse::
   :ref: stglib.core.cmd.lisstcdf2nc_parser
   :prog: runlisstcdf2nc.py
