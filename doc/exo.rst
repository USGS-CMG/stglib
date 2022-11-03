YSI EXO
*******

Data will generally be processed using a series of run scripts. The first script for each instrument type
depends on two :doc:`configuration files </config>`.

Instrument data to raw .cdf
===========================

Convert from exported .csv file from KOR software to a raw netCDF file with .cdf extension using ``runexocsv2cdf.py``.

runexocsv2cdf.py
----------------

.. argparse::
   :ref: stglib.core.cmd.exocsv2cdf_parser
   :prog: runexocsv2cdf.py

Raw .cdf to clean .nc
=====================

Convert the raw .cdf data into an CF-compliant netCDF file with .nc extension using ``runexocdf2nc.py``, optionally including :doc:`atmospheric correction </atmos>` of the pressure data.

runexocdf2nc.py
---------------

.. argparse::
   :ref: stglib.core.cmd.exocdf2nc_parser
   :prog: runexocdf2nc.py
