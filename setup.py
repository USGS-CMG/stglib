import pathlib

from setuptools import find_packages, setup

import versioneer

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / "README.md").read_text(encoding="utf-8")

# used https://github.com/pypa/sampleproject/blob/master/setup.py

setup(
    name="stglib",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description=("Process data from a variety of oceanographic instrumentation"),
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="U.S. Geological Survey",
    author_email="dnowacki@usgs.gov",
    url="https://github.com/USGS-CMG/stglib",
    license="Public domain",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    # install_requires=['numpy', 'netCDF4', 'xarray'],
    python_requires=">=3.8",
    packages=find_packages(exclude=["doc", "tests"]),
    scripts=[
        "scripts/runaqdhdr2cdf.py",
        "scripts/runaqdcdf2nc.py",
        "scripts/runhwlbcsv2cdf.py",
        "scripts/runhwlbcdf2nc.py",
        "scripts/runiqmat2cdf.py",
        "scripts/runiqcdf2nc.py",
        "scripts/runrskrsk2cdf.py",
        "scripts/runrskcsv2cdf.py",
        "scripts/runrskcdf2nc.py",
        "scripts/runrsknc2waves.py",
        "scripts/runrsknc2diwasp.py",
        "scripts/runecocsv2cdf.py",
        "scripts/runecocdf2nc.py",
        "scripts/runexocsv2cdf.py",
        "scripts/runexocdf2nc.py",
        "scripts/runtrollcsv2cdf.py",
        "scripts/runtrollcdf2nc.py",
        "scripts/runwvswad2cdf.py",
        "scripts/runwvscdf2nc.py",
        "scripts/runwvsnc2diwasp.py",
        "scripts/runwvsnc2waves.py",
        "scripts/aqdturnaround.py",
        "scripts/exoturnaround.py",
        "scripts/runwxtcsv2cdf.py",
        "scripts/runwxtcdf2nc.py",
        "scripts/runeofelog2cdf.py",
        "scripts/runeofecdf2nc.py",
        "scripts/runvecdat2cdf.py",
        "scripts/runveccdf2nc.py",
        "scripts/runsigmat2cdf.py",
        "scripts/runsigraw2cdf.py",
        "scripts/runsigcdf2nc.py",
    ],
    include_package_data=True,
)
