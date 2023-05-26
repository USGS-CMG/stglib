import os
import subprocess
import sysconfig
from pathlib import Path

import pytest

scripts = Path(sysconfig.get_path("scripts"))


def exo_raw(glob_att, config_yaml):
    result = subprocess.run(
        [scripts / "runexocsv2cdf.py", glob_att, config_yaml],
        capture_output=True,
        cwd="stglib/tests/data",
    )
    assert "Finished writing data" in result.stdout.decode("utf8")


def exo_nc(nc_file):
    result = subprocess.run(
        [scripts / "runexocdf2nc.py", nc_file],
        capture_output=True,
        cwd="stglib/tests/data",
    )
    assert "Done writing netCDF file" in result.stdout.decode("utf8")


def aqd_raw(glob_att, config_yaml):
    result = subprocess.run(
        [scripts / "runaqdhdr2cdf.py", glob_att, config_yaml],
        capture_output=True,
        cwd="stglib/tests/data",
    )
    assert "Finished writing data" in result.stdout.decode("utf8")


def aqd_nc(nc_file):
    result = subprocess.run(
        [scripts / "runaqdcdf2nc.py", nc_file],
        capture_output=True,
        cwd="stglib/tests/data",
    )
    assert "Done writing netCDF file" in result.stdout.decode("utf8")


def aqdhr_raw(glob_att, config_yaml):
    result = subprocess.run(
        [scripts / "runaqdhrhdr2cdf.py", glob_att, config_yaml],
        capture_output=True,
        cwd="stglib/tests/data",
    )
    assert "Finished writing data" in result.stdout.decode("utf8")


def aqdhr_nc(nc_file):
    result = subprocess.run(
        [scripts / "runaqdhrcdf2nc.py", nc_file],
        capture_output=True,
        cwd="stglib/tests/data",
    )
    assert "Done writing netCDF file" in result.stdout.decode("utf8")


def test_exo():
    exo_raw("glob_attbel5C.txt", "config_bel5C.yaml")
    exo_nc("bel53Cexo-raw.cdf")
    exo_raw("glob_att1119a.txt", "1119Aexo_config.yaml")
    exo_nc("1119Aexo-raw.cdf")
    exo_raw("glob_att1151b.txt", "1151Bexo_config.yaml")
    exo_nc("1151Bexo-raw.cdf")


def test_aqd():
    aqd_raw("glob_att1118a_b.txt", "aqd1118A_config.yaml")
    aqd_nc("1118ABaqd-raw.cdf")
    aqd_raw("glob_att1121a_msl_aqd.txt", "aqd1121A_config.yaml")
    aqd_nc("11211Aaqd-raw.cdf")


def test_aqdhr():
    aqdhr_raw("gatts_CHC14TDH.txt", "config_CHC14TDH.yaml")
    aqdhr_nc("CHC14TDH-raw.cdf")


def aqdturnaround(basefile):
    result = subprocess.run(
        ["python", "../../../scripts/aqdturnaround.py", basefile],
        capture_output=True,
        cwd="stglib/tests/data",
    )
    assert "Finished creating turnaround plots" in result.stdout.decode("utf8")


def test_aqdturnaround():
    # ensure plots are created for data collected in BEAM coordinates
    aqdturnaround("1121AQ04")
    # and XYZ coordinates
    aqdturnaround("NBMCCE02")
    # and ENU coordinates
    aqdturnaround("BEL503")


def wxt_raw(glob_att, config_yaml):
    result = subprocess.run(
        [scripts / "runwxtcsv2cdf.py", glob_att, config_yaml],
        capture_output=True,
        cwd="stglib/tests/data",
    )
    assert "Finished writing data" in result.stdout.decode("utf8")


def wxt_nc(nc_file):
    result = subprocess.run(
        [scripts / "runwxtcdf2nc.py", nc_file],
        capture_output=True,
        cwd="stglib/tests/data",
    )
    assert "Done writing netCDF file" in result.stdout.decode("utf8")


def test_wxt():
    wxt_raw("glob_att1149.txt", "wxt1149_config.yaml")
    wxt_nc("1149wxt-raw.cdf")


def iq_raw(glob_att, config_yaml):
    result = subprocess.run(
        [scripts / "runiqmat2cdf.py", glob_att, config_yaml],
        capture_output=True,
        cwd="stglib/tests/data",
    )
    assert "Finished writing data" in result.stdout.decode("utf8")


def iq_nc(nc_file):
    result = subprocess.run(
        [scripts / "runiqcdf2nc.py", nc_file],
        capture_output=True,
        cwd="stglib/tests/data",
    )
    assert "Done writing netCDF file" in result.stdout.decode("utf8")


def test_iq():
    iq_raw("glob_att1097C.txt", "config_1097C.yaml")
    iq_nc("10971Ciq-raw.cdf")


def eco_raw(glob_att, config_yaml):
    result = subprocess.run(
        [scripts / "runecocsv2cdf.py", glob_att, config_yaml],
        capture_output=True,
        cwd="stglib/tests/data",
    )
    assert "Finished writing data" in result.stdout.decode("utf8")


def eco_nc(nc_file):
    result = subprocess.run(
        [scripts / "runecocdf2nc.py", nc_file],
        capture_output=True,
        cwd="stglib/tests/data",
    )
    assert "Done writing netCDF file" in result.stdout.decode("utf8")


def test_eco():
    eco_raw("glob_att1103D.txt", "11032Decn_config.yaml")
    eco_nc("11032Decn-raw.cdf")


def rbr_raw(glob_att, config_yaml):
    result = subprocess.run(
        [scripts / "runrskcsv2cdf.py", glob_att, config_yaml],
        capture_output=True,
        cwd="stglib/tests/data",
    )
    assert "Finished writing data" in result.stdout.decode("utf8")


def rbr_nc(nc_file):
    result = subprocess.run(
        [scripts / "runrskcdf2nc.py", nc_file],
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
    rbr_raw("gatts_CSF20SC2.txt", "csf20sc201_config.yaml")
    rbr_nc("CSF20SC201pt-raw.cdf")


def eofe_raw(glob_att, config_yaml):
    result = subprocess.run(
        [scripts / "runeofelog2cdf.py", glob_att, config_yaml],
        capture_output=True,
        cwd="stglib/tests/data",
    )
    assert "Finished writing data" in result.stdout.decode("utf8")


def eofe_nc(nc_file):
    result = subprocess.run(
        [scripts / "runeofecdf2nc.py", nc_file],
        capture_output=True,
        cwd="stglib/tests/data",
    )
    assert "Done writing netCDF file" in result.stdout.decode("utf8")


def test_eofe():
    eofe_raw("glob_att1123A_msl.txt", "1123Aea_example_config.yaml")
    eofe_nc("11231Aea_example-raw.cdf")
    eofe_raw("glob_att1137.txt", "1137aa_config.yaml")
    eofe_nc("11373aa-raw.cdf")


def sig_mat(glob_att, config_yaml):
    result = subprocess.run(
        [scripts / "runsigmat2cdf.py", glob_att, config_yaml],
        capture_output=True,
        cwd="stglib/tests/data",
    )
    assert "Finished writing data" in result.stdout.decode("utf8")


def sig_nc(nc_file):
    result = subprocess.run(
        [scripts / "runsigcdf2nc.py", nc_file],
        capture_output=True,
        cwd="stglib/tests/data",
    )
    assert "Done writing netCDF file" in result.stdout.decode("utf8")


# @pytest.mark.skip(reason="works locally but not on github built-in checks")
def test_sig():
    sig_mat("glob_att1126_sig1.txt", "sig1126_config.yaml")
    print(os.listdir())
    sig_nc("11261sig_burst-raw.cdf")
    sig_mat("glob_att1126_sig2.txt", "sig11262_config.yaml")
    sig_nc("11262sig_burst-raw.cdf")
    sig_nc("11262sig_echo1-raw.cdf")


def ensure_cf(script, glob_att, config_yaml):
    result = subprocess.run(
        [scripts / script, glob_att, config_yaml],
        capture_output=True,
        cwd="stglib/tests/data",
    )
    assert "ValueError: Non-CF Conventions are not supported." in result.stderr.decode(
        "utf8"
    )
    assert result.returncode != 0


def test_ensure_cf():
    # ensure scripts fail if non-CF Conventions are specified
    ensure_cf(
        "runeofelog2cdf.py",
        "glob_att1123A_msl_EPIC.txt",
        "1123Aea_example_config.yaml",
    )
