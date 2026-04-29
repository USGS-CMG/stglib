from setuptools import setup

import versioneer

setup(
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    # install_requires=['numpy', 'netCDF4', 'xarray'],
    # packages=find_packages(exclude=["doc", "tests"]),
    scripts=[
        "scripts/aqdturnaround.py",
        "scripts/exoturnaround.py",
        "scripts/rbrturnaround.py",
    ],
    # include_package_data=True,
)
