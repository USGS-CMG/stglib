from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# used https://github.com/pypa/sampleproject/blob/master/setup.py

setup(name='stglib',
      version='0.1dev',
      description='Process oceanographic data',
      long_description=long_description,
      long_description_content_type='text/markdown',
      author='Dan Nowacki',
      author_email='dnowacki@usgs.gov',
      url='https://github.com/dnowacki-usgs/stglib',
      license='Public domain',
      classifiers=['Development Status :: 3 - Alpha',
                   'Intended Audience :: Science/Research',
                   'Programming Language :: Python',
                   'Programming Language :: Python :: 3',
                   'Programming Language :: Python :: 3.5',
                   'Programming Language :: Python :: 3.6',
                   'Programming Language :: Python :: 3.7'],
      # install_requires=['numpy', 'netCDF4', 'xarray'],
      python_requires='>=3.5',
      packages=find_packages(exclude=['doc', 'tests']),
      scripts=['scripts/runaqdhdr2cdf.py',
               'scripts/runaqdcdf2nc.py',
               'scripts/runhwlbcsv2cdf.py',
               'scripts/runhwlbcdf2nc.py',
               'scripts/runiqmat2cdf.py',
               'scripts/runiqcdf2nc.py',
               'scripts/runrskrsk2cdf.py',
               'scripts/runrskcdf2nc.py',
               'scripts/runrsknc2waves.py',
               'scripts/runrsknc2diwasp.py',
               'scripts/runecocsv2cdf.py',
               'scripts/runecocdf2nc.py',
               'scripts/runexocsv2cdf.py',
               'scripts/runexocdf2nc.py',
               'scripts/runwvswad2cdf.py',
               'scripts/runwvscdf2nc.py',
               'scripts/runwvsnc2diwasp.py',
               'scripts/runwvsnc2waves.py',
               'scripts/aqdturnaround.py',
               'scripts/exoturnaround.py'
               ],
      include_package_data=True
     )
