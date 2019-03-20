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
- ``<VAR>_min``: fill values less than this minimum valid value. Values outside this range will become ``_FillValue``. Substitute your variable for ``<VAR>``, e.g. ``fDOMQSU_min``.
- ``<VAR>_max``: fill values more than this maximum valid value.
- ``<VAR>_min_diff``: fill values where data decreases by more than this number of units in a single time step. Should be a negative number.
- ``<VAR>_max_diff``: fill values where data increases by more than this number of units in a single time step.
- ``<VAR>_med_diff``: fill values where difference between a 5-point (default) median filter and original values is greater than this number.
- ``<VAR>_med_diff_pct``: fill values where percent difference between a 5-point (default) median filter and original values is greater than this number.
- ``<VAR>_bad_ens``: specify bad ensembles that should be set to ``_FillValue``.
- ``trim_by_salinity``: if ``'true'``, use salinity (``S_41``) as a master variable. Wherever salinity is ``_FillValue``, all other variables will be filled as well. Useful for when the instrument comes out of the water.

.. literalinclude:: ../examples/exo_config.yaml
   :language: yaml
   :linenos:

NTU
---

NTU-specific options include:

- All the _min, _max, etc. options available to the EXO.
- ``Turb_std_max``: fill turbidity based on a maximum standard deviation value.
- ``spb``: samples per burst
- ``user_ntucal_coeffs``: polynomial coefficients, e.g., ``[9.078E-07, 5.883E-02, -2.899E+00]``.
