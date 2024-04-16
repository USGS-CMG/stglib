Vector
******

Instrument orientation
======================

Vector orientation refers to the direction of the z-axis, whether positive UP or DOWN.

.. note::
   For fixed stem vertical case systems, as well as cable probe vertical case systems, orientation = UP means the probe head is pointing DOWN (sample volume below probe). The I/O connector must be at the top, and probe cable at the bottom.

   For fixed stem vertical case systems, as well as cable probe vertical case systems, orientation = DOWN means the probe head is pointing UP (sample volume above probe). The I/O connector must be at the bottom, and probe cable at the top.

   For cable probe horizontal case systems, orientation = UP means the probe head is pointing UP (sample volume above probe) when the orientation marking on the case is UP.

   For cable probe horizontal case systems, orientation = DOWN means the probe head is pointing DOWN (sample volume above probe) when the orientation marking on the case is DOWN.

If your setup differs from any of the above, you may need to override the orientation with the ``orientation`` variable in the config file. This is a required variable and the user will be warned if it does not match what the instrument reports.

For more information, see "The Comprehensive Manual - Velocimeters" available from Nortek.

Instrument data to raw .cdf
===========================

First, export data from the Vector software to text format.

Convert from text to a raw netCDF file with ``.cdf`` extension using runvechdr2cdf.py. This script
depends on two arguments, the global attribute file and extra configuration information :doc:`configuration files </config>`.

runvecdat2cdf.py
----------------

.. argparse::
   :ref: stglib.core.cmd.vechdr2cdf_parser
   :prog: runvecdat2cdf.py

Raw .cdf to CF-compliant .nc
============================

Convert the raw .cdf data into an CF-compliant netCDF file with .nc extension, optionally including :doc:`atmospheric correction </atmos>` of the pressure data.

runveccdf2nc.py
---------------

.. argparse::
   :ref: stglib.core.cmd.veccdf2nc_parser
   :prog: runveccdf2nc.py

Compute wave statistics
=======================

Use stglib's built-in wave-statistics toolbox to compute a wave-statistics file.

**Experimental** PUV support is also present. Testing welcome!

runvecnc2waves.py
-----------------

.. argparse::
   :ref: stglib.core.cmd.vecnc2waves_parser
   :prog: runvecnc2waves.py
