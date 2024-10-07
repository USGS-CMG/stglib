SBE 26plus Seagauge (waves)
***************************
     
Data will generally be processed using a series of run scripts. The first script for each instrument type
depends on two :doc:`configuration files </config>`.

Instrument data to raw .cdf
===========================

Convert from exported .wb file to a raw netCDF file with .cdf extension using ``runsgwvswb2cdf.py``.

runsgwvswb2cdf.py
----------------

.. argparse::
   :ref: stglib.core.cmd.sgwvswb2cdf_parser
   :prog: runsgwvswb2cdf.py

Raw .cdf to CF-compliant .nc
============================

Convert the raw .cdf data into a CF-compliant netCDF file with .nc extension using ``runsgwvscdf2nc.py``.

runsgwvscdf2nc.py
---------------

.. argparse::
   :ref: stglib.core.cmd.sgwvscdf2nc_parser
   :prog: runsgwvscdf2nc.py
   
Compute wave statistics
=======================

Generate the waves statistics into a CF-compliant netCDF file with .nc extension using ``runsgwvsnc2waves.py``.

runsgwvsnc2waves.py
-----------------

.. argparse::
  :ref: stglib.core.cmd.sgwvsnc2waves_parser
  :prog: runsgwvsnc2waves.py
