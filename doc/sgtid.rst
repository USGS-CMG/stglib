Instrument data to raw .cdf
===========================

Convert from exported .tid file to a raw netCDF file with .cdf extension using ``runsgtid2cdf.py``.

runsgtid2cdf.py
----------------

.. argparse::
   :ref: stglib.core.cmd.sgtid2cdf_parser
   :prog: runsgtid2cdf.py

Raw .cdf to CF-compliant .nc
============================

Convert the raw .cdf data into a CF-compliant netCDF file with .nc extension using ``runsgcdf2nc.py``.

runsgcdf2nc.py
---------------

.. argparse::
   :ref: stglib.core.cmd.sgcdf2nc_parser
   :prog: runsgcdf2nc.py