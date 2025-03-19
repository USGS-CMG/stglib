Nortek Vector
*************

Instrument orientation
======================

Orientation in stglib refers to the probe head, whether pointing UP (sample volume above probe) or DOWN (sample volume below probe).

The orientation status code bit in the Vector data file refers to the direction of the z-axis, whether positive UP or down.
For more information, see "The Comprehensive Manual - Velocimeters" available from Nortek.

The ``orientation`` variable (indicating the probe head orientation) in the config file is required and the user will be warned if it does not match what the instrument reports.

Processing data
===============

First, export data from the Vector software to text format.

Then use :doc:`runots.py </runots>` to process using the two :doc:`configuration files </config>`.

Use stglib's built-in wave-statistics toolbox to compute a wave-statistics file (``runots.py vec nc2waves``).

**Experimental** PUV support is also present. Testing welcome!
