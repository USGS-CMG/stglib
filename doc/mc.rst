Seabird SBE 37 MicroCAT
**************

Data will generally be processed using a series of run scripts. The first script for each instrument type
depends on two :doc:`configuration files </config>`.

Instrument data to raw .cdf
===========================

Convert from exported .csv file to a raw netCDF file with .cdf extension using ``runmcasc2cdf.py``.

runmcasc2cdf.py
----------------

.. argparse::
   :ref: stglib.core.cmd.mcasc2cdf_parser
   :prog: runmcasc2cdf.py

Raw .cdf to CF-compliant .nc
============================

Convert the raw .cdf data into a CF-compliant netCDF file with .nc extension using ``runmccdf2nc.py``.

runmccdf2nc.py
---------------

.. argparse::
   :ref: stglib.core.cmd.mccdf2nc_parser
   :prog: runmccdf2nc.py
