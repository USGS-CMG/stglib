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
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    # install_requires=['numpy', 'netCDF4', 'xarray'],
    python_requires=">=3.9",
    packages=find_packages(exclude=["doc", "tests"]),
    scripts=[
        "scripts/aqdturnaround.py",
        "scripts/exoturnaround.py",
    ],
    entry_points={
        "console_scripts": [
            "runaqdhdr2cdf.py=stglib.core.runcmd:runaqdhdr2cdf",
            "runaqdcdf2nc.py=stglib.core.runcmd:runaqdcdf2nc",
            "runaqdhrhdr2cdf.py=stglib.core.runcmd:runaqdhrhdr2cdf",
            "runaqdhrcdf2nc.py=stglib.core.runcmd:runaqdhrcdf2nc",
            "runhobocsv2cdf.py=stglib.core.runcmd:runhobocsv2cdf",
            "runhobocdf2nc.py=stglib.core.runcmd:runhobocdf2nc",
            "runiqmat2cdf.py=stglib.core.runcmd:runiqmat2cdf",
            "runiqcdf2nc.py=stglib.core.runcmd:runiqcdf2nc",
            "runrskrsk2cdf.py=stglib.core.runcmd:runrskrsk2cdf",
            "runrskcsv2cdf.py=stglib.core.runcmd:runrskcsv2cdf",
            "runrskcdf2nc.py=stglib.core.runcmd:runrskcdf2nc",
            "runrsknc2waves.py=stglib.core.runcmd:runrsknc2waves",
            "runrsknc2diwasp.py=stglib.core.runcmd:runrsknc2diwasp",
            "runecocsv2cdf.py=stglib.core.runcmd:runecocsv2cdf",
            "runecocdf2nc.py=stglib.core.runcmd:runecocdf2nc",
            "runexocsv2cdf.py=stglib.core.runcmd:runexocsv2cdf",
            "runexocdf2nc.py=stglib.core.runcmd:runexocdf2nc",
            "runtrollcsv2cdf.py=stglib.core.runcmd:runtrollcsv2cdf",
            "runtrollcdf2nc.py=stglib.core.runcmd:runtrollcdf2nc",
            "runwvswad2cdf.py=stglib.core.runcmd:runwvswad2cdf",
            "runwvscdf2nc.py=stglib.core.runcmd:runwvscdf2nc",
            "runwvsnc2diwasp.py=stglib.core.runcmd:runwvsnc2diwasp",
            "runwvsnc2waves.py=stglib.core.runcmd:runwvsnc2waves",
            "runwxtcsv2cdf.py=stglib.core.runcmd:runwxtcsv2cdf",
            "runwxtcdf2nc.py=stglib.core.runcmd:runwxtcdf2nc",
            "runeofelog2cdf.py=stglib.core.runcmd:runeofelog2cdf",
            "runeofecdf2nc.py=stglib.core.runcmd:runeofecdf2nc",
            "runvecdat2cdf.py=stglib.core.runcmd:runvecdat2cdf",
            "runveccdf2nc.py=stglib.core.runcmd:runveccdf2nc",
            "runvecnc2waves.py=stglib.core.runcmd:runvecnc2waves",
            "runsigmat2cdf.py=stglib.core.runcmd:runsigmat2cdf",
            "runsigraw2cdf.py=stglib.core.runcmd:runsigraw2cdf",
            "runsigcdf2nc.py=stglib.core.runcmd:runsigcdf2nc",
            "runlisstcsv2cdf.py=stglib.core.runcmd:runlisstcsv2cdf",
            "runlisstcdf2nc.py=stglib.core.runcmd:runlisstcdf2nc",
            "runtcmcsv2cdf.py=stglib.core.runcmd:runtcmcsv2cdf",
            "runtcmcdf2nc.py=stglib.core.runcmd:runtcmcdf2nc",
            "runmcasc2cdf.py=stglib.core.runcmd:runmcasc2cdf",
            "runmccdf2nc.py=stglib.core.runcmd:runmccdf2nc",
        ],
    },
    include_package_data=True,
)
