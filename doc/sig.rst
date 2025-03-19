Nortek Signature
****************

**NOTE: this code works with up- or down-looking Signature data collected in 'beam', 'XYZ' or 'earth'.
It supports data types Burst, IBurst, Echo1, BurstHR, IBurstHR, Average. It also supports Altimeter and AHRS (Advanced Heading Reference System) data that are included when present with supported data types. Presently it does not yet support data type BottomTrack, but it can be added in the future as needed.**


First, export data from the Signature Deployment software to Matlab format with coordinate transformations.

Then use :doc:`runots.py </runots>` to process using the two :doc:`configuration files </config>`.

Use stglib's built-in wave-statistics toolbox and/or pyDIWASP to compute a wave-statistics file (``runots.py sig nc2waves`` , ``runots.py sig nc2diwasp``).
