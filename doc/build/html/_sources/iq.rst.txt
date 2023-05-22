SonTek IQ
*********

For data exploration purposes, stglib supports reading the ``.mat`` file exported from the SonTek IQ software into an xarray ``Dataset`` using :py:meth:`~stglib.iq.read_iq`.

Data will generally be processed using a series of run scripts. The first script for each instrument type
depends on two :doc:`configuration files </config>`.

IQ-exported Matlab file to raw .cdf
===================================

First, export .mat files from the instrument software.

Convert from exported .mat to a raw netCDF file with .cdf extension using ``runiqmat2cdf.py``.

runiqmat2cdf.py
----------------

.. argparse::
   :ref: stglib.core.cmd.iqmat2cdf_parser
   :prog: runiqmat2cdf.py

Raw .cdf to CF-compliant .nc
============================

Convert the raw .cdf data into an CF-compliant netCDF file with .nc extension using ``runiqcdf2nc.py``.

runiqcdf2nc.py
---------------

.. argparse::
   :ref: stglib.core.cmd.iqcdf2nc_parser
   :prog: runiqcdf2nc.py
