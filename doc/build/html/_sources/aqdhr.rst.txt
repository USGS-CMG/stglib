Aquadopp HR
***********

**Aquadopp HR mode data has only seen limited testing. Contributions welcome.**

Data are processed using a series of run scripts that use command line arguments.

Instrument data to raw .cdf
===========================

First, in AquaPro HR, export data to text files using the default options.

Convert from text to a raw netCDF file with ``.cdf`` extension using runaqdhrhdr2cdf.py. This script
depends on two arguments, the global attribute file and extra configuration information :doc:`configuration files </config>`.

runaqdhrhdr2cdf.py
------------------

.. argparse::
   :ref: stglib.core.cmd.aqdhdr2cdf_parser
   :prog: runaqdhrhdr2cdf.py

Raw .cdf to CF-compliant .nc
============================

Convert the raw .cdf data into an CF-compliant netCDF file with .nc extension, optionally including :doc:`atmospheric correction </atmos>` of the pressure data.  Correcting pressure for atmospheric is a side-bar task- use the .ipynb examples to see what to do.

runaqdhrcdf2nc.py
-----------------

.. argparse::
   :ref: stglib.core.cmd.aqdcdf2nc_parser
   :prog: runaqdhrcdf2nc.py
