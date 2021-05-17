import subprocess
import pytest


def exo_raw(glob_att, config_yaml):
    result = subprocess.run(
        ["runexocsv2cdf.py", glob_att, config_yaml],
        capture_output=True,
        cwd="stglib/tests/data",
    )
    assert "Finished writing data" in result.stdout.decode("utf8")


def exo_nc(nc_file):
    result = subprocess.run(
        ["runexocdf2nc.py", nc_file], capture_output=True, cwd="stglib/tests/data",
    )
    assert "Done writing netCDF file" in result.stdout.decode("utf8")


def test_exo():
    exo_raw("glob_attbel5C.txt", "config_bel5C.yaml")
    exo_nc("bel53Cexo-raw.cdf")
    exo_raw("glob_att1119a.txt", "1119Aexo_config.yaml")
    exo_nc("1119Aexo-raw.cdf")
