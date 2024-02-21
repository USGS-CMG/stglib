Aquadopp (waves)
****************

Data will generally be processed using a series of run scripts that use command line arguments.  For AQD waves it's a 3 step process.

Instrument data to raw .cdf
===========================

First, export data to text format using the AquaPro software and default options.

Convert from text to a raw netCDF file with ``.cdf`` extension using runwvswad2cdf.py. This script
depends on two arguments, the global attribute file and extra configuration information :doc:`configuration files </config>`.

runwvswad2cdf.py
----------------

.. argparse::
   :ref: stglib.core.cmd.wvswad2cdf_parser
   :prog: runwvswad2cdf.py


Raw .cdf to CF-compliant .nc
============================

Convert the raw .cdf data into an CF-compliant netCDF file with .nc extension, optionally including :doc:`atmospheric correction </atmos>` of the pressure data.

runwvscdf2nc.py
---------------

.. argparse::
  :ref: stglib.core.cmd.wvscdf2nc_parser
  :prog: runwvscdf2nc.py

Compute wave statistics
=======================

Using DIWASP: This is run separately in Matlab using the DIWASP toolbox (docs TODO).

Use stglib's built-in wave-statistics toolbox to compute a wave-statistics file.

runwvsnc2waves.py
-----------------

.. argparse::
  :ref: stglib.core.cmd.wvsnc2waves_parser
  :prog: runwvsnc2waves.py
