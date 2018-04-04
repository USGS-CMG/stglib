Configuration files
*******************

There are two required configuration files for processing data. Contents of both files will be included as attributes in both the xarray Dataset and the netCDF files.

Global attributes configuration file
====================================

This file describes attributes that apply to the mooring, and uses a peculiar formatting as shown in the example below.

.. literalinclude:: examples/glob_att1076a.txt
   :linenos:

Instrument-specific configuration file
======================================

This file is instrument-specific and is YAML formatted. An example is given below.

.. literalinclude:: examples/config.yaml
   :language: yaml
   :linenos:
