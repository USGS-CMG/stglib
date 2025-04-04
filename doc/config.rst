Configuration files
*******************

There are two required configuration files for processing data: the global attributes file, which describes attributes that apply to the mooring, and the instrument configuration file, which describes attributes that apply to an instrument on a mooring. Contents of both files will be included as attributes in both the xarray Dataset and the netCDF files.

A note on time and time zones
=============================

Time is always in Coordinated Universal Time (UTC).

Transitioning from EPIC to CF Conventions
=========================================

Historically, data have been released according to NOAA PMEL/EPIC conventions. Today, `CF Conventions <http://cfconventions.org>`_ are used much more frequently, and stglib supports only CF Conventions. Specifying conventions is done via ``Conventions`` keyword in either the global attributes file or the instrument configuration file.

Setting CF in global attributes
-------------------------------

::

  Conventions; CF-1.9

Setting CF in the instrument configuration file
-----------------------------------------------

::

  Conventions: 'CF-1.9'

Specifying CF-1.9 or a later release of the standard will enable straight-to-CF processing.

Global attributes configuration file
====================================

This file describes attributes that apply to the mooring, and uses a peculiar formatting as shown in the example below.

.. literalinclude:: ../examples/glob_att1076a.txt
   :linenos:

Instrument configuration file
=============================

This file is instrument-specific and is YAML formatted. A few examples are given below.

.. note::
   Although YAML supports boolean values, netCDF does not support them as attributes. Because stglib saves the values specified in the instrument configuration file as netCDF attributes, you must enclose values potentially interpreted as boolean (such as true or false) in quotation marks in the YAML file.

Options common to most (all?) instrument config files:

- ``Conventions``: version of the CF Conventions, ``'CF-1.9'`` presently
- ``basefile``: the input filename without extension
- ``filename``: output filename, to which ``-raw.cdf``, ``-a.nc``, etc. will be appended
- ``ClockError``: number, in seconds, negative is slow. Applies a simple offset for times. Useful if the instrument was deployed in an incorrect time zone.
- ``ClockDrift``: number, in seconds, negative is slow. Linearly interpolates times for when the instrument clock has drifted.
- ``initial_instrument_height``: elevation of instrument in meters
- ``initial_instrument_height_note``
- ``P_1ac_note``: a note on the atmospheric pressure source used
- ``zeroed_pressure``: a note detailing whether the pressure sensor was zeroed before deployment, and other pertinent details such as date and time of zeroing.
- ``good_dates``: a list of dates to clip data by instead of the default ``Deployment_date`` and ``Recovery_date``. Example: ``good_dates: ['2021-01-22 18:32', '2021-04-13 19:27'] # first burst looked suspect``. Multiple date ranges can also be used. Example: ``good_dates: ['2021-01-22 18:32', '2021-02-28 23:59', '2021-04-01 00:00', '2021-04-13 19:27'] # the month of March was bad``
- ``good_ens``: a list of good indices (based on the raw file, zero-based) to clip the data by. Example: ``good_ens: [10, 500]``. To specify multiple good ranges, add additional pairs of indices: ``good_ens: [10, 500, 560, 600]`` will clip the data to samples 10-500 and 560-600 in the final file.
- ``vert_dim``: user specified coordinate variable for vertical dimension for data variables with non-singular vertical dimension (default = 'z')

Multiple instruments
--------------------

Options applicable to many instrument types include:

- ``<VAR>_bad_ens``: specify bad ensemble ranges (either index numbers or dates) that should be set to ``_FillValue``. If you want multiple ranges, you can do this with additional values in the array. For example, ``Turb_bad_ens: ['2017-09-30 21:15', '2017-10-02 09:30', '2017-10-12 20:45', '2017-10-16 00:30']``. This will set the ranges in late September and early October, and again in mid-October, to ``_FillValue``.
- ``<VAR>_bad_ens_indiv``: specify ensembles (either index numbers or dates) that should be set to ``_FillValue``. For example, ``Turb_bad_ens_indiv: ['2017-09-30 21:15', '2017-10-02 09:30', '2017-10-12 20:45', '2017-10-16 00:30']``. This will set these four individual timestamps to ``_FillValue``.
- ``<VAR>_min``: fill values less than this minimum valid value. Values outside this range will become ``_FillValue``. Substitute your variable for ``<VAR>``, e.g. ``fDOMQSU_min``.
- ``<VAR>_max``: fill values more than this maximum valid value.
- ``<VAR>_min_diff``: fill values where data decreases by more than this number of units in a single time step. This will typically be a negative number.
- ``<VAR>_min_diff_pct``: fill values where data decreases by more than this percent in a single time step. This will typically be a negative number.
- ``<VAR>_max_diff``: fill values where data increases by more than this number of units in a single time step.
- ``<VAR>_max_diff_pct``: fill values where data increases by more than this percent in a single time step.
- ``<VAR>_med_diff``: fill values where difference between a 5-point (default) median filter and original values is greater than this number.
- ``<VAR>_med_diff_pct``: fill values where percent difference between a 5-point (default) median filter and original values is greater than this number.
- ``<VAR>_max_blip``: fill short-lived maximum "blips", values that increase greater than this number and then immediately decrease at the next time step.
- ``<VAR>_max_blip_pct``: fill short-lived maximum "blips", values that increase more than this percent and then immediately decrease at the next time step.
- ``<VAR>_trim_fliers``: fill flier values, which are data points surrounded by filled data. Set to the maximum size of flier clumps to remove.
- ``<VAR>_warmup_samples``: fill these many samples at the beginning of each burst.
- ``<VAR>_mask``: a single variable or list of variables which should be used to fill the given variable. For example ``u_1205_mask: ["cor1_1285", "cor2_1286", "cor3_1287"]`` will set ``u_1205`` to ``_FillValue`` wherever the correlation variables are ``_FillValue``
- ``<VAR>_mask_expr``: trim values based on an expression containing another variable. For example, ``Turb_mask_expr: "P_1ac < 0.1"`` will fill all ``Turb`` data where ``P_1ac`` is less than 0.1. Currently supported operators are ['>', '<', '>=', '<=', '==', '!=']. Usage is currently limited to simple expressions with the masking variable on the left-hand side.
- ``drop_vars``: a list of variables to be removed from the final file. For example, ``drop_vars: ['nLF_Cond_µS_per_cm', 'Wiper_Position_volt', 'Cable_Pwr_V']``.

Options for signal filtering:

- ``<VAR>_lowpass_filt``: apply butterworth lowpass filter with specified cutoff period in seconds.
- ``<VAR>_highpass_filt``: apply butterworth highpass filter with specified cutoff period in seconds.
- ``<VAR>_bandpass_filt``: apply butterworth bandpass filter with specified cutoff period in seconds as two element list [cut_long, cut_short].
- ``<VAR>_med_filt``: apply n point median filter, where n is specified value (must be an odd number).
- ``filter_order``: specify order of butterworth filter (default = 4 if not specified).


Options for wave processing using pyDIWASP:

- ``diwasp_method``: estimator method used by pyDIWASP (options (available now): 'IMLM' (default) or 'DFTM')
- ``diwasp_nfft``: length of FFTs used to calculate spectra (default = 256)
- ``diwasp_nsegs``: number of segments to use to window input data for spectral analysis (default = 16)
- ``diwasp_dres``: number (integer) of directions used for directional spectra (default = 180)
- ``diwasp_dunit``: specify directional units for directional wave spectra (options: 'naut' (default), 'cart', 'rad')
- ``diwasp_xdir``: The compass direction of the x-axis from which directions are measured (default = 90)
- ``diwasp_iter``: iteration limit (default = 50)
- ``diwasp_smooth``: smooth directional spectra using pyDIWASP smoothing option (default = 'ON')
- ``diwasp_ibin``: velocity bin number (0 = 1st bin) to use for directional wave processing (default = 0)
- ``diwasp_nsamps``: user specified number of samples to use in processing for each wave burst (optional)
- ``diwasp_pow2``: if set to 'true' use next lowest power of 2 of samples for processing each wave burst (default = 'false')

Refer to DIWASP original documentation for addition information: “DIWASP, a directional wave spectra toolbox for MATLAB®: User Manual. Research Report WP-1601-DJ (V1.4), Centre for Water Research, University of Western Australia.”

Aquadopp
--------

Aquadopp-specific options include:

- ``trim_method``: can be ``'water level'``, ``'water level sl'``, ``'bin range'``, ``None``, or ``'none'``. Or just omit the option entirely if you don't want to use it.
- ``<VAR>_trim_single_bins``: trim data where only a single bin of data (after trimming via ``trim_method``) remains. Set this value to ``true`` to enable.
- ``<VAR>_maxabs_diff_2d``: trim values in a 2D DataArray when the absolute value of the increase is greater than a specified amount
- ``AnalogInput1_<ATTR>`` or ``AnalogInput2_<ATTR>``: if ``<ATTR>`` is "height", "standard_name", "long_name", "units", "institution", "comment", "source", or "references", this will create the appropriate attribute for the given variable.

For Aquadopp waves:

- ``puv``: set to ``true`` to compute PUV wave statistics. **(EXPERIMENTAL)**

.. literalinclude:: ../examples/aqd_config.yaml
   :language: yaml
   :linenos:

Signature
---------

Signature-specific options include (see Aquadopp for others):

- ``outdir``: output directory (make sure it exists) to write individual ``cdf`` files before being compiled into a single ``cdf`` file per data type
- ``orientation``: can be ``UP`` or ``DOWN`` use this to identify orientation of profiler
- ``chunks``: list of key, value pairs for user specified chunking of data (e.g. ['time', 256000, 'bindist', 64])
- ``wave_interval``: interval in seconds for calculating wave bursts from continuous data
- ``wave_start_time``: start datetime for first wave burst (e.g. "2021-03-10 16:00:00")
- ``wp_min``, ``wp_max``: min/max allowable wave period, in seconds
- ``wh_min``, ``wh_max``: min/max allowable wave height, in meters
- ``wp_ratio``: maximum allowable ratio between peak period (``wp_peak``) and mean period (``wp_4060``).
- ``diwasp``: processing type for pyDIWASP wave processing; options available are: 'suv' and 'puv' for directional waves and 'pres' and 'elev' for non-directional waves
- ``puv``: if 'true' and ``nc2waves`` processing is called directional wave processing using stglib ``puv_quick_vectorized`` method is run in addition to standard non-directional stglib wave processing

.. literalinclude:: ../examples/aqd_config.yaml
   :language: yaml
   :linenos:

RBR instruments
---------------

Options specific to RBR instruments exported from the Ruskin software include:

- ``basefile``: the input filename without extension or data type. For example, if your exported text files are named ``055170_20190219_1547_burst.txt``, ``055170_20190219_1547_data.txt``, etc., ``basefile`` will be ``055170_20190219_1547``.
- ``filtered_wl``: "true" to turn on filtered water level variable (4th order lowpass butterworth filter with 6 min cutoff)
- ``wp_min``, ``wp_max``: min/max allowable wave period, in seconds
- ``wh_min``, ``wh_max``: min/max allowable wave height, in meters
- ``wp_ratio``: maximum allowable ratio between peak period (``wp_peak``) and mean period (``wp_4060``).
- ``<VAR>_min``: fill values less than this minimum valid value. Values outside this range will become ``_FillValue``. Substitute your variable for ``<VAR>``, e.g. ``P_1ac_min``. Only works for ``P_1`` and ``P_1ac``. Useful for trimming by minimum pressure for instruments that go dry on some tidal cycles. Any data within the burst less than the threshold will result in the full burst being filled.

.. literalinclude:: ../examples/dw_config.yaml
   :language: yaml
   :linenos:

When an RBR instrument is used in ``CONTINUOUS`` mode as a profiling instrument (e.g., twisting the endcap to start/stop a profile), include the following line in your configuration file:

- ``featureType: 'profile'``: this `CF-compliant <https://cfconventions.org/cf-conventions/cf-conventions.html#profile-data>`_ ``featureType`` instructs stglib to process these data as a profile dataset.
- ``latitude: [36.959, 41.533, 27.764]``, ``longitude: [-122.056, -70.651, -82.638]``: these values can each be specified as a YAML list of latitudes and longitudes, each element in the lists corresponding to a profile.
- ``split_profiles``: when set to `True`, split a multi-profile dataset into individual netCDF files for each profile

EXO
---

EXO-specific options include:

- ``skiprows``: number of lines to skip in the CSV before the real data begins

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

EofE ECHOLOGGER
---------------
- All the _min, _max, _bad_ens, etc. options available to the EXO.
- ``instrument_type``: types "ea" and "aa" are supported.
- ``orientation``: orientation of transducers types 'DOWN' or 'UP' are supported.
- ``average_salinity``: average salinity value (PSU) for the water mass for the deployment site and time period.
- ``average_salinity_note``: source of average salinity value.

Sequoia Scientific LISST
------------------------

- ``operating_mode``: set to ``burst`` if instrument was deployed in burst mode
- ``basefile_lop``: filename, without extension, of .lop metadata file (if present)

Sontek IQ
---------

- All the _min, _max, _bad_ens, etc. options available to the EXO.
- ``orientation``: can be ``UP`` or ``DOWN`` use this to identify orientation of profiler
- ``positive_direction``: direction (degrees) of positive flow indicated by the X arrow on top of instrument (optional, recommended)
- ``flood_direction``: direction (degrees) of flood current in channel, may be opposite of positive flow direction depending on field set up (optional, recommended)
- ``channel_cross_section_note``: note specifying starting bank (left or right) for RTK transect across the channel and when the transect measurements were collected (optional, recommended)

Onset Hobo
----------

- All the _min, _max, _bad_ens, etc. options available to the EXO.
- ``filtered_wl``: "true" to turn on filtered water level variable (4th order lowpass butterworth filter with 6 min cutoff)
- ``instrument_type``: can be ``hwl`` (water level), ``hwlb`` (water level as barometer), ``hdo`` (dissolved oxygen) or ``hcnd`` (conductivity) use these based on parameter measured by hobo logger
- ``skipfooter``: number of lines to skip in the CSV file at the end of the file
- ``ncols``: number of columns of data to read, starting at first
- ``names``: option for user specified column names (only recommended when code will not read names using automated/default method)

Lowell TCM Hobo
---------------

- All the _min, _max, _bad_ens, etc. options available to the EXO.
- ``skipfooter``: number of lines to skip in the CSV file at the end of the file
- ``ncols``: number of columns of data to read, starting at first
- ``names``: option for user specified column names (only recommended when code will not read names using automated/default method)

Vector
------
- ``pressure_sensor_height`` and ``velocity_sample_volume_height`` to specify the elevations of these two sensors.
- ``puv``: set to ``true`` to compute PUV wave statistics. **(EXPERIMENTAL)**
- ``orientation``: ``UP`` means probe head is pointing up (sample volume above probe head). ``DOWN`` means probe head is pointing down (sample volume below probe head).
- Many of the Aquadopp options apply to the Vector.

SBE 37 MicroCAT
---------------
- All the _min, _max, _bad_ens, etc. options available to the EXO
- ``skiprows``: number of lines to skip in the ASC before the real data begins

SBE 26plus Seagauge
-------------------
For Seagauge tides:
- ``file_type``: set to ``.tid`` or ``.wb`` to specify raw data file type
- ``calculated_tide_interval``: enter the desired tide interval when using the .wb file
- ``calculated_tide_interval_units``: tide interval units
- ``calculated_tide_duration``: enter the desired tide duration when using the .wb file
- ``calculated_tide_duration_units``: tide duration units
- All the _min, _max, _bad_ens, etc. options available to the EXO.

For Seagauge waves:
- ``calculated_wave_interval``: enter the desired wave interval
- ``calculated_wave_interval_units``: wave interval units
- ``wp_min``, ``wp_max``: min/max allowable wave period, in seconds
- ``wh_min``, ``wh_max``: min/max allowable wave height, in meters
- ``wp_ratio``: maximum allowable ratio between peak period (``wp_peak``) and mean period (``wp_4060``)

TruBlue
-------
- All the _min, _max, _bad_ens, etc. options available to the EXO
- ``skiprows``: number of header lines to skip in the txt file before the real data begins
- ``filtered_wl``: "true" to turn on filtered water level variable (4th order lowpass butterworth filter with 6 min cutoff)
- ``wave_interval``: interval in seconds for calculating wave bursts from continuous pressure data
- ``wp_min``, ``wp_max``: min/max allowable wave period, in seconds
- ``wh_min``, ``wh_max``: min/max allowable wave height, in meters
- ``wp_ratio``: maximum allowable ratio between peak period (``wp_peak``) and mean period (``wp_4060``).

AQUAscat1000R
-------------
- ``outdir``: path to desired folder for burst .cdf files (converted from burst .mat files in mat2cdf)
- ``basefile``: path to folder containing burst .mat files generated from Aquatec's ReadAquascat1000.m
- ``P_1_offset``: offset between 0 and abs pressure before deploying
- ``P_1_scale``: scale factor to apply to raw pressure data, likely 2. Plot raw data to asses appropriate scale factor.
- ``Tx_offset``: offset between ABSS temperature and actual temperature. Plot raw data to asses appropriate scale factor. Check with another instrument.
- ``Tx_scale``: scale factor to apply to raw temperature data. Likely none needed. Plot raw data to asses appropriate scale factor.
- ``Bat_offset``: offset between ABSS raw battery voltage and actual battery voltage. Likely none needed.
- ``Bat_scale``: scale factor to apply to raw battery voltage data, likely 2. Plot raw data to asses appropriate scale factor.
- ``orientation``: orientation of transducer(s)
- ``initial_instrument_height``: height of acoustic transducer
- ``pressure_sensor_height``: height of pressure port on canister, likely different that transducer height

Geolux Wave Radar
-----------------
- All the _min, _max, _bad_ens, etc. options available to the EXO
- ``filtered_wl``: "true" to turn on filtered water level variable (4th order lowpass butterworth filter with 6 min cutoff)
- ``wave_interval``: interval in seconds for calculating wave bursts from continuous sea-surface elevation data
- ``wave_duration``: duration in seconds for calculating wave statistics
- ``wp_min``, ``wp_max``: min/max allowable wave period, in seconds
- ``wh_min``, ``wh_max``: min/max allowable wave height, in meters
- ``wp_ratio``: maximum allowable ratio between peak period (``wp_peak``) and mean period (``wp_4060``).
- ``wavedat_tolerance``: tolerance in seconds to fill gaps in wave data to be used for calculating wave statistics (default = '20 s').
- ``wlfilt_tolerance``: tolerance in seconds to fill gaps in water level data to be used calculating fileterd water level ``water_level_filt`` (default = '60 s').

Teledyne RDI instruments
------------------------
- Mostly follows the options available to the Aquadopp
- Does not yet support waves
