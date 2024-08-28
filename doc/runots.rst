Processing data with runots
***************************

Data are processed using the umbrella run script ``runots``. With this script the user specifies the instrument and step to be applied to the data.

For most instruments, this is a two-step process:

1. instrument data to raw CDF
2. raw CDF to clean .nc

For some instruments (usually those measuring waves) there is a third step:

3. clean .nc to wave-statistics .nc

The user will call ``runots`` from the command line, with arguments as follows:

runots
------

.. argparse::
   :ref: stglib.core.cmd.runots_parser
   :prog: runots
