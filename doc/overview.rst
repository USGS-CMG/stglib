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
* Use :doc:`runots.py </runots>` to process the data.
* Most instruments require two or more steps, run in series, for full processing.
* When external data is required (e.g., for atmospheric compensation), these are provided to :doc:`runots.py </runots>`. In the case of atmospheric compensation, Jupyter notebooks are available to help create an appropriate atmospheric pressure record.
