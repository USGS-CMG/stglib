Processing RBR d|wave data
**************************

Data will generally be processed using a series of run scripts. The first script for each instrument type
depends on two :doc:`configuration files </config>`.

Instrument data to raw .cdf
===========================

Convert from .rsk binary to a raw netCDF file with .cdf extension using ``runrskrsk2cdf.py``.

runrskrsk2cdf.py
----------------

.. argparse::
   :ref: stglib.core.cmd.rskrsk2cdf_parser
   :prog: runrskrsk2cdf.py

Raw .cdf to clean .nc
=====================

Convert the raw .cdf data into an EPIC-compliant netCDF file with .nc extension using ``runrskcdf2nc.py``, optionally including :doc:`atmospheric correction </atmos>` of the pressure data.

runrskcdf2nc.py
---------------

.. argparse::
   :ref: stglib.core.cmd.rskcdf2nc_parser
   :prog: runrskcdf2nc.py


DIWASP
======

Run DIWASP (within MATLAB) to produce wave statistics (see ``scripts/rundiwasp.m`` for an example run script). DIWASP must be run within MATLAB.

Clean .nc to DIWASP .nc
=======================

Incorporate the DIWASP statistics into an EPIC-compliant netCDF file with .nc extension using ``runrsknc2diwasp.py``.

runrsknc2diwasp.py
------------------

.. argparse::
  :ref: stglib.core.cmd.rsknc2diwasp_parser
  :prog: runrsknc2diwasp.py
