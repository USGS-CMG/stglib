Configuration files
*******************

There are two required configuration files for processing data. Contents of both files will be included as attributes in both the xarray Dataset and the netCDF files.

Global attributes configuration file
====================================

This file describes attributes that apply to the mooring, and uses a peculiar formatting as shown in the example below.

.. literalinclude:: ../examples/glob_att1076a.txt
   :linenos:

Instrument configuration file
=============================

This file is instrument-specific and is YAML formatted. A few examples are given below.

Options common to most (all?) instrument config files:

- ``basefile``: the input filename without extension
- ``filename``: output filename, to which ``-raw.cdf``, ``-a.nc``, etc. will be appended
- ``LatLonDatum``: will likely be ``'NAD83'``. TODO: should this be in glob_att instead?
- ``ClockError``: number, in seconds, negative is slow.
- ``initial_instrument_height``: elevation of instrument in meters
- ``initial_instrument_height_note``
- ``P_1ac_note``: a note on the atmospheric pressure source used

Aquadopp
--------

Aquadopp-specific options include:

- ``head_rotation``: probably will be ``'horizontal'``
- ``zeroed_pressure``
- ``cutoff_ampl``: will probably always be ``0``
- ``trim_method``: can be ``'water level'``, ``'water level sl'``, ``'bin range'``, ``None``, or ``'none'``. Or just omit the option entirely if you don't want to use it.

.. literalinclude:: ../examples/aqd_config.yaml
   :language: yaml
   :linenos:

d|wave
------

d|wave-specific options include:

- ``wp_min``, ``wp_max``: min/max allowable wave period, in seconds
- ``wh_min``, ``wh_max``: min/max allowable wave height, in meters
- ``wp_ratio``: maximum allowable ratio between peak period (``wp_peak``) and mean period (``wp_4060``).

.. literalinclude:: ../examples/dw_config.yaml
   :language: yaml
   :linenos:

EXO
---

EXO-specific options include:

- ``skiprows``: number of lines to skip in the CSV before the real data begins
- ``fDOMRFU_max_diff``: maximum point-to-point difference between consecutive values
- ``C_51_min_diff``: minimum point-to-point difference between consecutive values
- ``fDOMQSU_max_diff``: each variable has a ``_max_diff`` and ``_min_diff`` option
- ``SpC_48_min_diff``
- ``S_41_min_diff``
- ``Turb_max_diff``

.. literalinclude:: ../examples/exo_config.yaml
   :language: yaml
   :linenos:
