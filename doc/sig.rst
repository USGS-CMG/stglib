Signature
*********

**NOTE: this code works with up- or down-looking Signature data collected in 'beam', 'XYZ' or 'earth'.
It supports data types Burst, IBurst, Echo1, BurstHR, and IBurstHR. It also supports Altimeter and AHRS (Advanced Heading Reference System) data that are included when present with supported data types. Presently it does not yet support data types Average and BottomTrack, but they will be added in the future as needed.**


First, export data from the Signature Deployment software to Matlab format with coordinate transformations.

Then use :doc:`runots.py </runots>` to process using the two :doc:`configuration files </config>`.
