Nortek instruments
******************

Aquadopp (currents)
===================

First, in AquaPro, export data to text files using the default options.

Then use :doc:`runots.py </runots>` to process using the two :doc:`configuration files </config>`.

Aquadopp (waves)
================

First, export data to text format using the AquaPro software and default options.

Then use :doc:`runots.py </runots>` to process using the two :doc:`configuration files </config>`.

Use stglib's built-in wave-statistics toolbox to compute a wave-statistics file (``runots.py wvs nc2waves``).

**Experimental** PUV support is also present. Testing welcome!

Aquadopp HR
===========

**Aquadopp HR mode data has only seen limited testing. Contributions welcome.**

First, in AquaPro HR, export data to text files using the default options.

Then use :doc:`runots.py </runots>` to process using the two :doc:`configuration files </config>`.

Signature
=========

**NOTE: this code works with up- or down-looking Signature data collected in 'beam', 'XYZ' or 'earth'.
It supports data types Burst, IBurst, Echo1, BurstHR, IBurstHR, Average. It also supports Altimeter and AHRS (Advanced Heading Reference System) data that are included when present with supported data types. Presently it does not yet support data type BottomTrack, but it can be added in the future as needed.**


First, export data from the Signature Deployment software to Matlab format with coordinate transformations.

Then use :doc:`runots.py </runots>` to process using the two :doc:`configuration files </config>`.

In the cdf2nc, nc2waves, and nc2diwasp steps, you can optionally input a .nc file that contains salinity and water temperature to be used to improve the estimates of water_level and water_depth variables by including [--salwtemp SALWTEMP]. Note: SALWTEMP nc file needs to include time period of the deployment being processed.

Use stglib's built-in wave-statistics toolbox and/or pyDIWASP to compute a wave-statistics file (``runots.py sig nc2waves`` , ``runots.py sig nc2diwasp``).

Nortek Vector
=============

Instrument orientation
----------------------

Orientation in stglib refers to the probe head, whether pointing UP (sample volume above probe) or DOWN (sample volume below probe).

The orientation status code bit in the Vector data file refers to the direction of the z-axis, whether positive UP or down.
For more information, see "The Comprehensive Manual - Velocimeters" available from Nortek.

The ``orientation`` variable (indicating the probe head orientation) in the config file is required and the user will be warned if it does not match what the instrument reports.

Processing data
---------------

First, export data from the Vector software to text format.

Then use :doc:`runots.py </runots>` to process using the two :doc:`configuration files </config>`.

Use stglib's built-in wave-statistics toolbox to compute a wave-statistics file (``runots.py vec nc2waves``).

**Experimental** PUV support is also present. Testing welcome!
