Configuration files
*******************

There are two required configuration files for processing data. Contents of both files will be included as attributes in both the xarray Dataset and the netCDF files.

Transitioning from EPIC to CF Conventions
=========================================

Historically, data have been released according to NOAA PMEL/EPIC conventions. Today, `CF Conventions <http://cfconventions.org>`_ are used much more frequently, and moving forward most data sets will be released following CF. This is supported within stglib by use of the ``Conventions`` keyword in the instrument configuration file. Setting this to:

::

  Conventions: 'CF-1.6'

(or a later release of the standard) will enable straight-to-CF processing.

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
- ``zeroed_pressure``: a note detailing whether the pressure sensor was zeroed before deployment, and other pertinent details such as date and time of zeroing.
- ``good_dates``: a list of dates to clip data by instead of the default ``Deployment_date`` and ``Recovery_date``. Example: ``good_dates: ['2021-01-22 18:32', '2021-04-13 19:27'] # first burst looked suspect``
- ``good_ens``: a list of good indices (based on the raw file, zero-based) to clip the data by. Example: ``good_ens: [10, 500]``. To specify multiple good ranges, add additional pairs of indices: ``good_ens: [10, 500, 560, 600]`` will clip the data to samples 10-500 and 560-600 in the final file.

Aquadopp
--------

Aquadopp-specific options include:

- ``head_rotation``: probably will be ``'horizontal'``
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
- ``<VAR>_min``: fill values less than this minimum valid value. Values outside this range will become ``_FillValue``. Substitute your variable for ``<VAR>``, e.g. ``P_1ac_min``. Only works for ``P_1`` and ``P_1ac``. Useful for trimming by minimum pressure for instruments that go dry on some tidal cycles. Any data within the burst less than the threshold will result in the full burst being filled.

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
- ``<VAR>_bad_ens``: specify bad ensemble ranges (either index numbers or dates) that should be set to ``_FillValue``. If you want multiple ranges, you can do this with additional values in the array. For example, ``Turb_bad_ens: ['2017-09-30 21:15', '2017-10-02 09:30', '2017-10-12 20:45', '2017-10-16 00:30']``. This will set the ranges in late September and early October, and again in mid-October, to ``_FillValue``.
- ``trim_by_salinity``: if ``'true'``, use salinity (``S_41``) as a master variable. Wherever salinity is ``_FillValue``, all other variables will be filled as well. Useful for when the instrument comes out of the water.
- ``drop_vars``: a list of variables to be removed from the final file. For example, ``drop_vars: ['nLF_Cond_µS_per_cm', 'Wiper_Position_volt', 'Cable_Pwr_V']``.

Note that negative numeric values in the YAML config file must be treated with care so as not to be interpreted as strings. If you want the minimum value to be, say, -0.2 units for a particular parameter, you must write this as ``-0.2`` and not ``-.2`` in the config file. The latter format will be interpreted as a string and will cause an error.

.. literalinclude:: ../examples/exo_config.yaml
   :language: yaml
   :linenos:

WET Labs ECO NTU
----------------

NTU-specific options include:

- All the _min, _max, _bad_ens, etc. options available to the EXO.
- ``Turb_std_max``: fill turbidity based on a maximum standard deviation value.
- ``spb``: samples per burst
- ``user_ntucal_coeffs``: polynomial coefficients, e.g., ``[9.078E-07, 5.883E-02, -2.899E+00]``.

Vaisala WXT536
--------------

WXT-specific options include:

- ``RTK_elevation_NAVD88``: RTK elevation of the sensor referenced to NAVD88 in meters.
- ``dir_offset``: a direction offset in degrees from magnetic north to be applied if the sensor was not pointing toward magnetic north.
- ``dir_offset_note``: a note about the direction offset being used.