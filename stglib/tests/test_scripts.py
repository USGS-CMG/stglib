import subprocess
import pytest


def exo_raw(glob_att, config_yaml):
    result = subprocess.run(
        ["python", "../../../scripts/runexocsv2cdf.py", glob_att, config_yaml],
        capture_output=True,
        cwd="stglib/tests/data",
    )
    assert "Finished writing data" in result.stdout.decode("utf8")


def exo_nc(nc_file):
    result = subprocess.run(
        ["python", "../../../scripts/runexocdf2nc.py", nc_file],
        capture_output=True,
        cwd="stglib/tests/data",
    )
    assert "Done writing netCDF file" in result.stdout.decode("utf8")


def aqd_raw(glob_att, config_yaml):
    result = subprocess.run(
        ["python", "../../../scripts/runaqdhdr2cdf.py", glob_att, config_yaml],
        capture_output=True,
        cwd="stglib/tests/data",
    )
    assert "Finished writing data" in result.stdout.decode("utf8")


def aqd_nc(nc_file):
    result = subprocess.run(
        ["python", "../../../scripts/runaqdcdf2nc.py", nc_file],
        capture_output=True,
        cwd="stglib/tests/data",
    )
    assert "Done writing netCDF file" in result.stdout.decode("utf8")

def test_exo():
    exo_raw("glob_attbel5C.txt", "config_bel5C.yaml")
    exo_nc("bel53Cexo-raw.cdf")
    exo_raw("glob_att1119a.txt", "1119Aexo_config.yaml")
    exo_nc("1119Aexo-raw.cdf")

def test_aqd():
    aqd_raw("glob_att1118a_b.txt", "aqd1118A_config.yaml")
    aqd_nc("1118ABaqd-raw.cdf")
    
def wxt_raw(glob_att, config_yaml):
    result = subprocess.run(
        ["python", "../../../scripts/runwxtcsv2cdf.py", glob_att, config_yaml],
        capture_output=True,
        cwd="stglib/tests/data",
    )
    assert "Finished writing data" in result.stdout.decode("utf8")

def wxt_nc(nc_file):
    result = subprocess.run(
        ["python", "../../../scripts/runwxtcdf2nc.py", nc_file],
        capture_output=True,
        cwd="stglib/tests/data",
    )
    assert "Done writing netCDF file" in result.stdout.decode("utf8")
    
def test_wxt():
    wxt_raw("glob_att1149.txt", "wxt1149_config.yaml")
    wxt_nc("1149wxt-raw.cdf")