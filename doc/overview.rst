Processing overview
*******************

Data from many instrument types can be processed using stglib. The general paradigm of processing is that each instrument is associated with a mooring (which could be a bottom platform, a wire mooring, a tripod, or something else). A mooring may have more than one instrument.


The general procedure is as follows:

For each mooring
================

* Create a global attributes :doc:`configuration file </config>`.

For each instrument
===================

* Create an instrument-specific :doc:`configuration file </config>`.
* Use the appropriate run script(s) to process the data.
* Some instruments require the use of two or more run scripts, run in series, for full processing.
* When external data is required (e.g., for atmospheric compensation), these are provided to the run scripts. In the case of atmospheric compensation, Jupyter notebooks are available to help create an appropriate atmospheric pressure record.
