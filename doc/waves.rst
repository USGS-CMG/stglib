Processing waves
****************

stglib supports the computing of wave statistics, both using pressure as well as PUV for directional wave statistics.

Waves can be computed both with fixed frequency cutoffs, as well as dynamic cutoffs based on water depth.

Spectra tails beyond the cutoff frequency are applied following `Jones & Monismith (2007) <https://doi.org/10.4319/lom.2007.5.317>`_.

For more information, see the code in `stglib/core/waves.py`
