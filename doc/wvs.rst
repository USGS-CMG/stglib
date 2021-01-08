Processing Aquadopp (waves) data
***********************************

Data will generally be processed using a series of run scripts that use command line arguments.  For AQD waves it's a 3 step process.

Step 1 : Instrument data to raw .cdf
=====================================

Step 1 : Convert from text to a raw netCDF file with ``.cdf`` extension using runaqdhdr2cdf.py. This script
depends on two arguments, the global attribute file and extra configuration information :doc:`configuration files </config>`.

runwvswad2cdf.py
----------------

.. argparse::
   :ref: stglib.core.cmd.wvswad2cdf_parser
   :prog: runwvswad2cdf.py


Step 2 : Convert the raw .cdf to clean, EPIC format .nc using runaqdcdf2nc.py.
==============================================================================

Convert the raw .cdf data into an EPIC-compliant netCDF file with .nc extension, optionally including :doc:`atmospheric correction </atmos>` of the pressure data.  Correcting pressure for atmospheric is a side-bar task- use the .ipynb examples to see what to do.

runwvscdf2nc.py
---------------

.. argparse::
  :ref: stglib.core.cmd.wvscdf2nc_parser
  :prog: runwvscdf2nc.py

Step 3 : Compute waves statistics
=================================

Using DIWASP: This is run separately in Matlab using the DIWASP toolbox (docs TODO).

Use stglib's built-in waves statistics toolbox to compute a wave statistics file.

runwvsnc2waves.py
-----------------

.. argparse::
  :ref: stglib.core.cmd.wvsnc2waves_parser
  :prog: runwvsnc2waves.py
