Processing Aquadopp data
************************

Data will generally be processed using a series of run scripts. The first script for each instrument type
depends on two :doc:`configuration files </config>`.

Instrument data to raw .cdf
===========================

.. argparse::
   :ref: stglib.core.cmd.aqdhdr2cdf_parser
   :prog: runaqdhdr2cdf.py

Raw .cdf to clean .nc
=====================

.. argparse::
   :ref: stglib.core.cmd.aqdcdf2nc_parser
   :prog: runaqdcdf2nc.py
