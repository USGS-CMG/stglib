Processing waves
****************

stglib supports the computing of wave statistics, both using pressure as well as PUV for directional wave statistics.
Waves can be computed both with fixed frequency cutoffs, as well as dynamic cutoffs based on water depth.
Spectra tails beyond the cutoff frequency are applied following `Jones & Monismith (2007) <JM>`_.
For more information, see the code in ``stglib/core/waves.py``.

In brief, stglib's wave-statistics code (:py:func:`stglib.core.waves.make_waves_ds`) does the following:

#. Compute the pressure spectra using Welch's method in :func:`stglib.core.waves.pressure_spectra`.
#. Compute the transfer function from wavenumber, water depth, and sensor height in :func:`stglib.core.waves.transfer_function`.
#. Compute the surface-elevation spectra by dividing the pressure spectra by the transfer function in :meth:`stglib.core.waves.elevation_spectra`.
#. Define the cutoff following `Jones & Monismith (2007) <JM>`_ in :func:`stglib.core.waves.define_cutoff`.
#. Add a tail following `Jones & Monismith (2007) <JM>`_ in :func:`stglib.core.waves.make_tail`.
#. Compute the zeroth and second moments of the surface-elevation spectra using :func:`stglib.core.waves.make_moment`.
#. Compute significant wave height, mean period, and peak period using :func:`stglib.core.waves.make_Hs`, :func:`stglib.core.waves.make_Tm`, and :func:`stglib.core.waves.make_Tp`.

The above list is for information only. The user does not need to apply these steps manually; they are all called by the run script.


.. autosummary::
  :toctree: generated/

  stglib.core.waves.make_waves_ds
  stglib.core.waves.pressure_spectra
  stglib.core.waves.transfer_function
  stglib.core.waves.elevation_spectra
  stglib.core.waves.define_cutoff
  stglib.core.waves.make_tail
  stglib.core.waves.make_moment
  stglib.core.waves.make_Hs
  stglib.core.waves.make_Tm
  stglib.core.waves.make_Tp

.. _JM: https://doi.org/10.4319/lom.2007.5.317
