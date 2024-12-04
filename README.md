# stglib - Process data from a variety of oceanographic instrumentation

[doi:10.5066/P13IQYFW](https://doi.org/10.5066/P13IQYFW)

This package contains code to process data from a variety of oceanographic instrumentation, consistent with the procedures of the USGS [Coastal/Marine Hazards and Resources Program](https://marine.usgs.gov) (formerly Coastal and Marine Geology Program).

Currently, this package has at least partial support for:

- Nortek Aquadopp profilers, in mean-current wave-burst, and HR modes
- Nortek Vector velocimeters
- Nortek Signature profilers
- RBR pressure (including waves) and turbidity sensors
- YSI EXO2 water-quality sondes
- SonTek IQ flow monitors
- WET labs sensors, including ECO NTUSB and ECO PAR
- Onset HOBO pressure sensors
- Vaisala Weather Transmitter WXT sensors
- In-Situ Aqua TROLL sensors
- RD Instruments ADCPs
- Moving-boat ADCP data processed using [QRev](https://hydroacoustics.usgs.gov/movingboat/QRev.shtml), for use in index-velocity computation
- EofE ECHOLOGGER altimeters
- SBE 37 Microcat
- SBE 26plus Seagauge
- TruBlue pressure sensors

We have plans to support:

- RDI Sentinel V profilers

This package makes heavy use of [NumPy](http://www.numpy.org), [xarray](http://xarray.pydata.org/en/stable/), and [netCDF4](http://unidata.github.io/netcdf4-python/). It works on Python 3.9+.

[Read the documentation](http://stglib.readthedocs.io/).
