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


def iq_raw(glob_att, config_yaml):
    result = subprocess.run(
        ["python", "../../../scripts/runiqmat2cdf.py", glob_att, config_yaml],
        capture_output=True,
        cwd="stglib/tests/data",
    )
    assert "Finished writing data" in result.stdout.decode("utf8")


def iq_nc(nc_file):
    result = subprocess.run(
        ["python", "../../../scripts/runiqcdf2nc.py", nc_file],
        capture_output=True,
        cwd="stglib/tests/data",
    )
    assert "Done writing netCDF file" in result.stdout.decode("utf8")


def test_iq():
    iq_raw("glob_att1097C.txt", "config_1097C.yaml")
    iq_nc("10971Ciq-raw.cdf")


def eco_raw(glob_att, config_yaml):
    result = subprocess.run(
        ["python", "../../../scripts/runecocsv2cdf.py", glob_att, config_yaml],
        capture_output=True,
        cwd="stglib/tests/data",
    )
    assert "Finished writing data" in result.stdout.decode("utf8")


def eco_nc(nc_file):
    result = subprocess.run(
        ["python", "../../../scripts/runecocdf2nc.py", nc_file],
        capture_output=True,
        cwd="stglib/tests/data",
    )
    assert "Done writing netCDF file" in result.stdout.decode("utf8")


def test_eco():
    eco_raw("glob_att1103D.txt", "11032Decn_config.yaml")
    eco_nc("11032Decn-raw.cdf")


def rbr_raw(glob_att, config_yaml):
    result = subprocess.run(
        ["python", "../../../scripts/runrskcsv2cdf.py", glob_att, config_yaml],
        capture_output=True,
        cwd="stglib/tests/data",
    )
    assert "Finished writing data" in result.stdout.decode("utf8")


def rbr_nc(nc_file):
    result = subprocess.run(
        ["python", "../../../scripts/runrskcdf2nc.py", nc_file],
        capture_output=True,
        cwd="stglib/tests/data",
    )
    assert "Done writing netCDF file" in result.stdout.decode("utf8")


def test_rbr():
    import zipfile
    from os import path

    # zip files created on windows won't extract to directories on unix systems without the following line
    path.altsep = "\\"
    with zipfile.ZipFile("stglib/tests/data/051001_CSF20SC201.zip", "r") as zip_ref:
        zip_ref.extractall("stglib/tests/data/")
    rbr_raw("csf20sc201_globatt.txt", "csf20sc201_config.yaml")
    rbr_nc("CSF20SC201pt-raw.cdf")

def eofe_raw(glob_att, config_yaml):
    result = subprocess.run(
        ["python", "../../../scripts/runeofelog2cdf.py", glob_att, config_yaml],
        capture_output=True,
        cwd="stglib/tests/data",
    )
    assert "Finished writing data" in result.stdout.decode("utf8")

def eofe_nc(nc_file):
    result = subprocess.run(
        ["python", "../../../scripts/runeofecdf2nc.py", nc_file],
        capture_output=True,
        cwd="stglib/tests/data",
    )
    assert "Done writing netCDF file" in result.stdout.decode("utf8")


def test_eofe():
    eofe_raw("glob_att1123A_msl.txt", "1123Aea_example_config.yaml")
    eofe_nc("11231Aea_example-raw.cdf")

