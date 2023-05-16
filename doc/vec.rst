Vector
******

Data will generally be processed using a series of run scripts that use command line arguments.

Instrument data to raw .cdf
===========================

First, export data from the Vector software to text format.

Convert from text to a raw netCDF file with ``.cdf`` extension using runvechdr2cdf.py. This script
depends on two arguments, the global attribute file and extra configuration information :doc:`configuration files </config>`.

runvechdr2cdf.py
----------------

.. argparse::
   :ref: stglib.core.cmd.vechdr2cdf_parser
   :prog: runvecdat2cdf.py

Raw .cdf to CF-compliant .nc
============================

Convert the raw .cdf data into an CF-compliant netCDF file with .nc extension, optionally including :doc:`atmospheric correction </atmos>` of the pressure data.  Correcting pressure for atmospheric is a side-bar task- use the .ipynb examples to see what to do.

runveccdf2nc.py
---------------

.. argparse::
   :ref: stglib.core.cmd.veccdf2nc_parser
   :prog: runveccdf2nc.py
