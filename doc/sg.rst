Seabird SBE 26plus Seagauge
***************************
     
Data will generally be processed using a series of run scripts. The first script for each instrument type
depends on two :doc:`configuration files </config>`.

Instrument data to raw .cdf
===========================

Convert from exported .tide file to a raw netCDF file with .cdf extension using ``runsgtid2cdf.py``.

runsgtid2cdf.py
----------------

.. argparse::
   :ref: stglib.core.cmd.sgtid2cdf_parser
   :prog: runsgtid2cdf.py

Raw .cdf to CF-compliant .nc
============================

Convert the raw .cdf data into a CF-compliant netCDF file with .nc extension using ``runsgcdf2nc.py``.

runsgcdf2nc.py
---------------

.. argparse::
   :ref: stglib.core.cmd.sgcdf2nc_parser
   :prog: runsgcdf2nc.py
