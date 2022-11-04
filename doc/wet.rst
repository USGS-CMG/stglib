WET Labs ECO
************

Data will generally be processed using a series of run scripts. The first script for each instrument type
depends on two :doc:`configuration files </config>`.

Instrument data to raw .cdf
===========================

Convert from  tab-delimited log file to a raw netCDF file with .cdf extension using ``runecocsv2cdf.py``.

runecocsv2cdf.py
----------------

.. argparse::
   :ref: stglib.core.cmd.ecolog2cdf_parser
   :prog: runecocsv2cdf.py

Raw .cdf to clean .nc
=====================

Convert the raw .cdf data into a CF-compliant netCDF file with .nc extension using ``runecocdf2nc.py``.

runecocdf2nc.py
---------------

.. argparse::
   :ref: stglib.core.cmd.ecocdf2nc_parser
   :prog: runecocdf2nc.py
