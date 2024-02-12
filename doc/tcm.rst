Lowell Tilt Current Meter (TCM)
*******************************

Data will generally be processed using a series of run scripts. The first script for each instrument type
depends on two :doc:`configuration files </config>`.

Instrument data to raw .cdf
===========================

Convert from exported CSV file to a raw netCDF file with .cdf extension using ``runtcmcsv2cdf.py``.

runtcmcsv2cdf.py
-----------------

.. argparse::
   :ref: stglib.core.cmd.tcmcsv2cdf_parser
   :prog: runtcmcsv2cdf.py

Raw .cdf to CF-compliant .nc
============================

Convert the raw .cdf data into CF-compliant netCDF file with .nc extension using ``runtcmcdf2nc.py``.

runtcmcdf2nc.py
----------------

.. argparse::
   :ref: stglib.core.cmd.tcmcdf2nc_parser
   :prog: runtcmcdf2nc.py
