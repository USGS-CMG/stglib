import subprocess
import sysconfig
import zipfile
from os import path
from pathlib import Path

import pytest

# zip files created on windows won't extract to directories on unix systems without the following line
path.altsep = "\\"

scripts = Path(sysconfig.get_path("scripts"))
cwd = "stglib/tests/data"


def exo_raw(glob_att, config_yaml):
    result = subprocess.run(
        [scripts / "runexocsv2cdf.py", glob_att, config_yaml],
        capture_output=True,
        cwd=cwd,
    )
    assert "Finished writing data" in result.stdout.decode("utf8")


def exo_nc(nc_file, atmpres=None):
    if atmpres is not None:
        runlist = [scripts / "runexocdf2nc.py", nc_file, "--atmpres", atmpres]
    else:
        runlist = [scripts / "runexocdf2nc.py", nc_file]
    result = subprocess.run(
        runlist,
        capture_output=True,
        cwd=cwd,
    )
    assert "Done writing netCDF file" in result.stdout.decode("utf8")


def test_exo():
    exo_raw("glob_attbel5C.txt", "config_bel5C.yaml")
    exo_nc("bel53Cexo-raw.cdf")
    exo_raw("glob_att1119a.txt", "1119Aexo_config.yaml")
    exo_nc("1119Aexo-raw.cdf")
    exo_raw("glob_att1151b.txt", "1151Bexo_config.yaml")
    exo_nc("1151Bexo-raw.cdf")
    # test for atmospheric correction
    exo_raw("glob_attbel3C.txt", "config_bel3C.yaml")
    exo_nc("BEL19B3C03exo-raw.cdf", "atmpres_BEL19B3C03exo.cdf")


def aqd_raw(glob_att, config_yaml):
    result = subprocess.run(
        [scripts / "runaqdhdr2cdf.py", glob_att, config_yaml],
        capture_output=True,
        cwd=cwd,
    )
    assert "Finished writing data" in result.stdout.decode("utf8")


def aqd_nc(nc_file, atmpres=None):
    if atmpres is not None:
        runlist = [scripts / "runaqdcdf2nc.py", nc_file, "--atmpres", atmpres]
    else:
        runlist = [scripts / "runaqdcdf2nc.py", nc_file]
    result = subprocess.run(
        runlist,
        capture_output=True,
        cwd=cwd,
    )
    assert "Done writing netCDF file" in result.stdout.decode("utf8")


def aqdhr_raw(glob_att, config_yaml):
    result = subprocess.run(
        [scripts / "runaqdhrhdr2cdf.py", glob_att, config_yaml],
        capture_output=True,
        cwd=cwd,
    )
    assert "Finished writing data" in result.stdout.decode("utf8")


def aqdhr_nc(nc_file):
    result = subprocess.run(
        [scripts / "runaqdhrcdf2nc.py", nc_file],
        capture_output=True,
        cwd=cwd,
    )
    assert "Done writing netCDF file" in result.stdout.decode("utf8")


def aqd_wvs_raw(glob_att, config_yaml):
    result = subprocess.run(
        [scripts / "runwvswad2cdf.py", glob_att, config_yaml],
        capture_output=True,
        cwd=cwd,
    )
    assert "Finished writing data" in result.stdout.decode("utf8")


def aqd_wvs_nc(nc_file, atmpres=None):
    if atmpres is not None:
        runlist = [scripts / "runwvscdf2nc.py", nc_file, "--atmpres", atmpres]
    else:
        runlist = [scripts / "runwvscdf2nc.py", nc_file]
    result = subprocess.run(
        runlist,
        capture_output=True,
        cwd=cwd,
    )
    assert "Done writing netCDF file" in result.stdout.decode("utf8")


def aqd_wvs_wvs(nc_file):
    result = subprocess.run(
        [scripts / "runwvsnc2waves.py", nc_file],
        capture_output=True,
        cwd=cwd,
    )
    assert "Done creating" in result.stdout.decode("utf8")


def test_aqd():
    aqd_raw("glob_att1118a_b.txt", "aqd1118A_config.yaml")
    aqd_nc("1118ABaqd-raw.cdf")
    aqd_raw("glob_att1121a_msl_aqd.txt", "aqd1121A_config.yaml")
    aqd_nc("11211Aaqd-raw.cdf")
    # test for atmospheric correction
    aqd_raw("nbm22cce01_gatts.txt", "config_nbm22cce01.yaml")
    aqd_nc("NBM22CCEaqd-raw.cdf", "atmpres_NBM22CCE.cdf")


def test_aqdhr():
    aqdhr_raw("gatts_CHC14TDH.txt", "config_CHC14TDH.yaml")
    aqdhr_nc("CHC14TDH-raw.cdf")
    aqdhr_raw("glob_att1113_aqdHR_tst.txt", "aqdhr1113tst_config.yaml")
    aqdhr_nc("1113aqdHR-raw.cdf")


def test_aqd_wvs():
    aqd_wvs_raw("glob_attbel5C.txt", "BEL5C_wvsconfig.yaml")
    aqd_wvs_nc("BEL19B5C04aqdwv-raw.cdf", atmpres="atmpres-BEL5Cwvs.cdf")
    aqd_wvs_wvs("BEL19B5C04aqdwvb-cal.nc")


def aqdturnaround(basefile):
    result = subprocess.run(
        ["python", "../../../scripts/aqdturnaround.py", basefile],
        capture_output=True,
        cwd=cwd,
    )
    assert "Finished creating turnaround plots" in result.stdout.decode("utf8")


def test_aqdturnaround():
    # ensure plots are created for data collected in BEAM coordinates
    aqdturnaround("1121AQ04")
    # and XYZ coordinates
    aqdturnaround("NBMCCE02")
    # and ENU coordinates
    aqdturnaround("BEL503")


def vec_raw(glob_att, config_yaml):
    result = subprocess.run(
        [scripts / "runvecdat2cdf.py", glob_att, config_yaml],
        capture_output=True,
        cwd=cwd,
    )
    assert "Finished writing data" in result.stdout.decode("utf8")


def vec_nc(nc_file, atmpres=None):
    if atmpres is not None:
        runlist = [scripts / "runveccdf2nc.py", nc_file, "--atmpres", atmpres]
    else:
        runlist = [scripts / "runveccdf2nc.py", nc_file]
    result = subprocess.run(
        runlist,
        capture_output=True,
        cwd=cwd,
    )
    assert "Done writing netCDF file" in result.stdout.decode("utf8")


def vec_wvs(nc_file):
    result = subprocess.run(
        [scripts / "runvecnc2waves.py", nc_file],
        capture_output=True,
        cwd=cwd,
    )
    assert "Done writing netCDF file" in result.stdout.decode("utf8")


def test_vec():
    vec_raw("gatts_NBM22CSB.txt", "config_NBM22CSB.yaml")
    vec_nc("NBMCSBvec01-raw.cdf")
    vec_wvs("NBMCSBvec01b-cal.nc")


def wxt_raw(glob_att, config_yaml):
    result = subprocess.run(
        [scripts / "runwxtcsv2cdf.py", glob_att, config_yaml],
        capture_output=True,
        cwd=cwd,
    )
    assert "Finished writing data" in result.stdout.decode("utf8")


def wxt_nc(nc_file):
    result = subprocess.run(
        [scripts / "runwxtcdf2nc.py", nc_file],
        capture_output=True,
        cwd=cwd,
    )
    assert "Done writing netCDF file" in result.stdout.decode("utf8")


def test_wxt():
    wxt_raw("glob_att1149.txt", "wxt1149_config.yaml")
    wxt_nc("1149wxt-raw.cdf")


def iq_raw(glob_att, config_yaml):
    result = subprocess.run(
        [scripts / "runiqmat2cdf.py", glob_att, config_yaml],
        capture_output=True,
        cwd=cwd,
    )
    assert "Finished writing data" in result.stdout.decode("utf8")


def iq_nc(nc_file):
    result = subprocess.run(
        [scripts / "runiqcdf2nc.py", nc_file],
        capture_output=True,
        cwd=cwd,
    )
    assert "Done writing netCDF file" in result.stdout.decode("utf8")


def test_iq():
    iq_raw("glob_att1097C.txt", "config_1097C.yaml")
    iq_nc("10971Ciq-raw.cdf")


def eco_raw(glob_att, config_yaml):
    result = subprocess.run(
        [scripts / "runecocsv2cdf.py", glob_att, config_yaml],
        capture_output=True,
        cwd=cwd,
    )
    assert "Finished writing data" in result.stdout.decode("utf8")


def eco_nc(nc_file):
    result = subprocess.run(
        [scripts / "runecocdf2nc.py", nc_file],
        capture_output=True,
        cwd=cwd,
    )
    assert "Done writing netCDF file" in result.stdout.decode("utf8")


def test_eco():
    eco_raw("glob_att1103D.txt", "11032Decn_config.yaml")
    eco_nc("11032Decn-raw.cdf")


def rbr_raw(glob_att, config_yaml):
    result = subprocess.run(
        [scripts / "runrskcsv2cdf.py", glob_att, config_yaml],
        capture_output=True,
        cwd=cwd,
    )
    assert "Finished writing data" in result.stdout.decode("utf8")


def rbr_nc(nc_file, atmpres=None):
    if atmpres is not None:
        runlist = [scripts / "runrskcdf2nc.py", nc_file, "--atmpres", atmpres]
    else:
        runlist = [scripts / "runrskcdf2nc.py", nc_file]
    result = subprocess.run(
        runlist,
        capture_output=True,
        cwd=cwd,
    )
    assert "Done writing netCDF file" in result.stdout.decode("utf8")


def rbr_wvs(nc_file):
    result = subprocess.run(
        [scripts / "runrsknc2waves.py", nc_file],
        capture_output=True,
        cwd=cwd,
    )
    assert "Done writing netCDF file" in result.stdout.decode("utf8")


def test_rbr():
    with zipfile.ZipFile("stglib/tests/data/051001_CSF20SC201.zip", "r") as zip_ref:
        zip_ref.extractall("stglib/tests/data/")
    rbr_raw("gatts_CSF20SC2.txt", "csf20sc201_config.yaml")
    rbr_nc("CSF20SC201pt-raw.cdf")
    rbr_wvs("CSF20SC201ptb-cal.nc")
    with zipfile.ZipFile("stglib/tests/data/055109_20220808_1605.zip", "r") as zip_ref:
        zip_ref.extractall("stglib/tests/data/")
    rbr_raw("gatts_055109_20220808_1605.txt", "055109_20220808_1605_config.yaml")
    rbr_nc("11512Cdw-raw.cdf")
    rbr_wvs("11512Cdwcont-cal.nc")


def test_rbr_profile():
    # profiling CTD
    with zipfile.ZipFile("stglib/tests/data/205598_20220104_1336.zip", "r") as zip_ref:
        zip_ref.extractall("stglib/tests/data/")
    rbr_raw(
        "gatts_205598_20220104_1336_wtw21.txt", "config_205598_20220104_1336_wtw21.yaml"
    )
    rbr_nc("WTW21CTD-raw.cdf", "atmpres_205598_20220104_1336_wtw21.cdf")


def eofe_raw(glob_att, config_yaml):
    result = subprocess.run(
        [scripts / "runeofelog2cdf.py", glob_att, config_yaml],
        capture_output=True,
        cwd=cwd,
    )
    assert "Finished writing data" in result.stdout.decode("utf8")


def eofe_nc(nc_file):
    result = subprocess.run(
        [scripts / "runeofecdf2nc.py", nc_file],
        capture_output=True,
        cwd=cwd,
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
        cwd=cwd,
    )
    assert "Finished writing data" in result.stdout.decode("utf8")


def sig_nc(nc_file):
    result = subprocess.run(
        [scripts / "runsigcdf2nc.py", nc_file],
        capture_output=True,
        cwd=cwd,
    )
    assert "Done writing netCDF file" in result.stdout.decode("utf8")


@pytest.mark.skip(reason="works locally but not on github built-in checks")
def test_sig():
    sig_mat("glob_att1126_sig1.txt", "sig1126_config.yaml")
    sig_nc("11261sig_burst-raw.cdf")
    sig_mat("glob_att1126_sig2.txt", "sig11262_config.yaml")
    sig_nc("11262sig_burst-raw.cdf")
    sig_nc("11262sig_echo1-raw.cdf")
    sig_mat("gatts_MIA23SH2_cf_rev.txt", "sig_avg_config.yaml")
    sig_nc("MIAsig_avgd-raw.cdf")
    sig_nc("MIAsig_altavgd-raw.cdf")


def hobo_raw(glob_att, config_yaml):
    result = subprocess.run(
        [scripts / "runhobocsv2cdf.py", glob_att, config_yaml],
        capture_output=True,
        cwd=cwd,
    )
    assert "Finished writing data" in result.stdout.decode("utf8")


def hobo_nc(nc_file):
    result = subprocess.run(
        [scripts / "runhobocdf2nc.py", nc_file],
        capture_output=True,
        cwd=cwd,
    )
    assert "Done writing netCDF file" in result.stdout.decode("utf8")


def test_hobo():
    hobo_raw("glob_att1168_hobowl.txt", "1168hwl_config.yaml")
    hobo_nc("11681hwl-raw.cdf")
    hobo_raw("glob_att1171_hobowl_baro.txt", "1171hwl_baro_config.yaml")
    hobo_nc("11711hwlb-raw.cdf")


def lisst_raw(glob_att, config_yaml):
    result = subprocess.run(
        [scripts / "runlisstcsv2cdf.py", glob_att, config_yaml],
        capture_output=True,
        cwd=cwd,
    )
    assert "Finished writing data" in result.stdout.decode("utf8")


def lisst_nc(nc_file):
    result = subprocess.run(
        [scripts / "runlisstcdf2nc.py", nc_file],
        capture_output=True,
        cwd=cwd,
    )
    assert "Done writing netCDF file" in result.stdout.decode("utf8")


def test_lisst():
    lisst_raw("gatts_lisst_L0221705.txt", "config_lisst_L0221705.yaml")
    lisst_nc("lisst_L0221705-raw.cdf")


def tcm_raw(glob_att, config_yaml):
    result = subprocess.run(
        [scripts / "runtcmcsv2cdf.py", glob_att, config_yaml],
        capture_output=True,
        cwd=cwd,
    )
    assert "Finished writing data" in result.stdout.decode("utf8")


def tcm_nc(nc_file):
    result = subprocess.run(
        [scripts / "runtcmcdf2nc.py", nc_file],
        capture_output=True,
        cwd=cwd,
    )
    assert "Done writing netCDF file" in result.stdout.decode("utf8")


def test_tcm():
    tcm_raw("glob_att1170_tcm.txt", "1170tcm_config.yaml")
    tcm_nc("11701tcm-raw.cdf")


def ensure_cf(script, glob_att, config_yaml):
    result = subprocess.run(
        [scripts / script, glob_att, config_yaml],
        capture_output=True,
        cwd=cwd,
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
