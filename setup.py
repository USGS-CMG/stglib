from setuptools import setup

setup(name='stglib',
      version='0.1dev',
      description='Process oceanographic data',
      author='Dan Nowacki',
      author_email='dnowacki@usgs.gov',
      url='https://github.com/dnowacki-usgs/stglib',
      license='Public domain',
      install_requires=['numpy', 'netCDF4', 'xarray'],
      packages=['stglib', 'stglib.aqd', 'stglib.rsk'],
      scripts=['scripts/runaqdhdr2cdf.py',
               'scripts/runaqdcdf2nc.py',
               'scripts/runrskrsk2cdf.py',
               'scripts/runrskcdf2nc.py',
               'scripts/runrsknc2diwasp.py',
               ],
      include_package_data=True
     )
