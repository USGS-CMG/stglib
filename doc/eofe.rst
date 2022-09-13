Processing EofE ECHOLOGGER data
*******************************

Data will generally be processed using a series of run scripts. The first script for each instrument type
depends on two :doc:`configuration files </config>`.

Instrument data to raw .cdf
===========================

Convert from  tab or space-delimited log file to a raw netCDF file with .cdf extension using ``runeofelog2cdf.py``.

runeofelog2cdf.py
-----------------

.. argparse::
   :ref: stglib.core.cmd.eofelog2cdf_parser
   :prog: runeofelog2cdf.py

Raw .cdf to clean .nc
=====================

Convert the raw .cdf data into CF-compliant netCDF file with .nc extension using ``runeofecdf2nc.py``.

runeofecdf2nc.py
----------------

.. argparse::
   :ref: stglib.core.cmd.eofecdf2nc_parser
   :prog: runeofecdf2nc.py
