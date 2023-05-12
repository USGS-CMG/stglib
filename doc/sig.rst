Signature
******

Data will generally be processed using a series of run scripts that use command line arguments.

**NOTE: this code works with up- or down-looking Signature data collected in 'beam', 'XYZ' or 'earth'.
It supports data types Burst, IBurst, Echo1, BurstHR, and IBurstHR. It also supports Altimeter and AHRS (Advanced Heading Reference System) data that are included when present with supported data types. Presently it does not yet support data types Average and BottomTrack, but they will be added in the future as needed.**

Instrument data to raw .cdf
===========================

First, export data from the Signature Deployment software to Matlab format with coordinate transformations.

Convert from mat files to a raw netCDF file (for each data type) with ``.cdf`` extension using runsigmat2cdf.py. This script
depends on two arguments, the global attribute file and extra configuration information :doc:`configuration files </config>`.

runsigmat2cdf.py
----------------

.. argparse::
   :ref: stglib.core.cmd.sigmat2cdf_parser
   :prog: runsigmat2cdf.py

Raw .cdf to CF-compliant .nc
============================

Convert the raw .cdf data (for each data type) into an CF-compliant netCDF file with .nc extension, optionally including :doc:`atmospheric correction </atmos>` of the pressure data.  Correcting pressure for atmospheric is a side-bar task- use the .ipynb examples to see what to do.

runsigcdf2nc.py
---------------

.. argparse::
   :ref: stglib.core.cmd.sigcdf2nc_parser
   :prog: runsigcdf2nc.py
