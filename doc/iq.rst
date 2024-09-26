SonTek IQ
*********

For data exploration purposes, stglib supports reading the ``.mat`` file exported from the SonTek IQ software into an xarray ``Dataset`` using :py:meth:`~stglib.iq.read_iq`.

Data will generally be processed using a series of run scripts. The first script for each instrument type
depends on two :doc:`configuration files </config>`.

IQ-exported Matlab file to raw .cdf
===================================

First, export .mat files from the instrument software. Then, use :doc:`runots </runots>` to process using the two :doc:`configuration files </config>`.
