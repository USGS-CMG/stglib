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


def run_script(script_name, *args):
    """Utility function to run a script with specified arguments and check output."""
    result = subprocess.run(
        [scripts / script_name] + list(args),
        capture_output=True,
        cwd=cwd,
    )
    print(result.stdout)
    output = result.stdout.decode("utf8")
    assert (
        "Finished writing data" in output
        or "Done writing netCDF file" in output
        or "Done creating" in output
        or "Finished creating" in output
        or "Done writing averaged netCDF file" in output
    )


def abss_raw(glob_att, config_yaml):
    run_script("runots.py", "abss", "mat2cdf", glob_att, config_yaml)


def abss_nc(nc_file, atmpres=None):
    if atmpres is not None:
        run_script("runots.py", "abss", "cdf2nc", nc_file, "--atmpres", atmpres)
    else:
        run_script("runots.py", "abss", "cdf2nc", nc_file)


def test_abss():
    abss_raw("glob_att1126_abs_test_msl.txt", "config_1126abs910_abs_test.yaml")
    abss_nc("11266abs910-131_test-raw.cdf")


def exo_raw(glob_att, config_yaml):
    run_script("runexocsv2cdf.py", glob_att, config_yaml)


def exo_nc(nc_file, atmpres=None):
    if atmpres is not None:
        run_script("runexocdf2nc.py", nc_file, "--atmpres", atmpres)
    else:
        run_script("runexocdf2nc.py", nc_file)


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
    run_script("runaqdhdr2cdf.py", glob_att, config_yaml)


def aqd_nc(nc_file, atmpres=None):
    if atmpres is not None:
        run_script("runaqdcdf2nc.py", nc_file, "--atmpres", atmpres)
    else:
        run_script("runaqdcdf2nc.py", nc_file)


def aqdhr_raw(glob_att, config_yaml):
    run_script("runaqdhrhdr2cdf.py", glob_att, config_yaml)


def aqdhr_nc(nc_file):
    run_script("runaqdhrcdf2nc.py", nc_file)


def aqd_wvs_raw(glob_att, config_yaml):
    run_script("runwvswad2cdf.py", glob_att, config_yaml)


def aqd_wvs_nc(nc_file, atmpres=None):
    if atmpres is not None:
        run_script("runwvscdf2nc.py", nc_file, "--atmpres", atmpres)
    else:
        run_script("runwvscdf2nc.py", nc_file)


def aqd_wvs_wvs(nc_file):
    run_script("runwvsnc2waves.py", nc_file)


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
    if not Path("stglib/tests/data/BEL503.hdr").is_file():
        with zipfile.ZipFile("stglib/tests/data/BEL503.zip", "r") as zip_ref:
            zip_ref.extractall("stglib/tests/data/")
    aqd_wvs_raw("glob_attbel5C.txt", "BEL5C_wvsconfig.yaml")
    aqd_wvs_nc("BEL19B5C04aqdwv-raw.cdf", atmpres="atmpres-BEL5Cwvs.cdf")
    aqd_wvs_wvs("BEL19B5C04aqdwvb-cal.nc")


def aqdturnaround(basefile):
    run_script("aqdturnaround.py", basefile)


def test_aqdturnaround():
    # ensure plots are created for data collected in BEAM coordinates
    aqdturnaround("1121AQ04")
    # and XYZ coordinates
    aqdturnaround("NBMCCE02")
    # and ENU coordinates
    if not Path("stglib/tests/data/BEL503.hdr").is_file():
        with zipfile.ZipFile("stglib/tests/data/BEL503.zip", "r") as zip_ref:
            zip_ref.extractall("stglib/tests/data/")
    aqdturnaround("BEL503")


def vec_raw(glob_att, config_yaml):
    run_script("runvecdat2cdf.py", glob_att, config_yaml)


def vec_nc(nc_file, atmpres=None):
    if atmpres is not None:
        run_script("runveccdf2nc.py", nc_file, "--atmpres", atmpres)
    else:
        run_script("runveccdf2nc.py", nc_file)


def vec_wvs(nc_file):
    run_script("runvecnc2waves.py", nc_file)


# @pytest.mark.skip(reason="works locally but not on GitLab CI")
def test_vec_burst():
    # burst mode vector, fails on gitlab CI
    with zipfile.ZipFile("stglib/tests/data/FISH0104green.zip", "r") as zip_ref:
        zip_ref.extractall("stglib/tests/data/")
    vec_raw("gatts_FISH0104green.txt", "FISH0104green.yaml")
    vec_nc("FISH0104green-raw.cdf")
    vec_wvs("FISH0104greenb.nc")


def test_vec_continuous():
    # continuous mode Vector
    with zipfile.ZipFile("stglib/tests/data/V1482304.zip", "r") as zip_ref:
        zip_ref.extractall("stglib/tests/data/")
    vec_raw("glob_att1126_msl.txt", "config_1126vec14823.yaml")
    vec_nc("1126vec14823-raw.cdf")


def wxt_raw(glob_att, config_yaml):
    run_script("runwxtcsv2cdf.py", glob_att, config_yaml)


def wxt_nc(nc_file):
    run_script("runwxtcdf2nc.py", nc_file)


def test_wxt():
    wxt_raw("glob_att1149.txt", "wxt1149_config.yaml")
    wxt_nc("1149wxt-raw.cdf")


def iq_raw(glob_att, config_yaml):
    run_script("runiqmat2cdf.py", glob_att, config_yaml)


def iq_nc(nc_file):
    run_script("runiqcdf2nc.py", nc_file)


def test_iq():
    iq_raw("glob_att1097C.txt", "config_1097C.yaml")
    iq_nc("10971Ciq-raw.cdf")


def eco_raw(glob_att, config_yaml):
    run_script("runecocsv2cdf.py", glob_att, config_yaml)


def eco_nc(nc_file):
    run_script("runecocdf2nc.py", nc_file)


def test_eco():
    eco_raw("glob_att1103D.txt", "11032Decn_config.yaml")
    eco_nc("11032Decn-raw.cdf")


def rbr_raw(glob_att, config_yaml):
    run_script("runrskcsv2cdf.py", glob_att, config_yaml)


def rbr_nc(nc_file, atmpres=None):
    if atmpres is not None:
        run_script("runrskcdf2nc.py", nc_file, "--atmpres", atmpres)
    else:
        run_script("runrskcdf2nc.py", nc_file)


def rbr_wvs(nc_file):
    run_script("runrsknc2waves.py", nc_file)


def rbr_diwasp(nc_file):
    run_script("runots.py", "rbr", "nc2diwasp", nc_file)


def test_rbr_tu():
    with zipfile.ZipFile("stglib/tests/data/054215_DMP23X3E01tu.zip", "r") as zip_ref:
        zip_ref.extractall("stglib/tests/data/")
    rbr_raw("gatts_DMP23X3E.txt", "DMP23X3E01tu_config.yaml")
    rbr_nc("DMP23X3E01tu-raw.cdf")


def test_rbr_wvs():
    with zipfile.ZipFile("stglib/tests/data/051001_CSF20SC201.zip", "r") as zip_ref:
        zip_ref.extractall("stglib/tests/data/")
    rbr_raw("gatts_CSF20SC2.txt", "csf20sc201_config.yaml")
    rbr_nc("CSF20SC201pt-raw.cdf")
    rbr_wvs("CSF20SC201ptb.nc")


def test_rbr():
    with zipfile.ZipFile("stglib/tests/data/055109_20220808_1605.zip", "r") as zip_ref:
        zip_ref.extractall("stglib/tests/data/")
    rbr_raw("gatts_055109_20220808_1605.txt", "055109_20220808_1605_config.yaml")
    rbr_nc("11512Cdw-raw.cdf")
    rbr_wvs("11512Cdwb.nc")
    # rbr_diwasp("11512Cdwcont-cal.nc") #need P_1c variable for this test to work


def test_rbr_profile():
    # profiling CTD
    with zipfile.ZipFile("stglib/tests/data/205598_20220104_1336.zip", "r") as zip_ref:
        zip_ref.extractall("stglib/tests/data/")
    rbr_raw(
        "gatts_205598_20220104_1336_wtw21.txt", "config_205598_20220104_1336_wtw21.yaml"
    )
    rbr_nc("WTW21CTD-raw.cdf", "atmpres_205598_20220104_1336_wtw21.cdf")


def eofe_raw(glob_att, config_yaml):
    run_script("runeofelog2cdf.py", glob_att, config_yaml)


def eofe_nc(nc_file):
    run_script("runeofecdf2nc.py", nc_file)


def test_eofe():
    eofe_raw("glob_att1123A_msl.txt", "1123Aea_example_config.yaml")
    eofe_nc("11231Aea_example-raw.cdf")
    eofe_raw("glob_att1137.txt", "1137aa_config.yaml")
    eofe_nc("11373aa-raw.cdf")


def sig_mat(glob_att, config_yaml):
    run_script("runots.py", "sig", "mat2cdf", glob_att, config_yaml)


def sig_nc(nc_file):
    run_script("runots.py", "sig", "cdf2nc", nc_file)


def sig_wvs(nc_file):
    run_script("runots.py", "sig", "nc2waves", nc_file)


def sig_diwasp(nc_file):
    run_script("runots.py", "sig", "nc2diwasp", nc_file)


@pytest.mark.skip(reason="works locally but not on github built-in checks")
def test_sig():
    sig_mat("glob_att1126_sig1.txt", "sig1126_config.yaml")
    sig_nc("11261sig_burst-raw.cdf")
    sig_mat("glob_att1126_sig2.txt", "sig11262_config.yaml")
    sig_nc("11262sig_echo1-raw.cdf")
    sig_nc("11262sig_burst-raw.cdf")
    sig_mat("gatts_MIA23SH2_cf_rev.txt", "sig_avg_config.yaml")
    sig_nc("MIAsig_avgd-raw.cdf")
    sig_nc("MIAsig_altavgd-raw.cdf")


@pytest.mark.skip(reason="works locally but not on github built-in checks")
def test_sig_wvs():
    sig_wvs("11261sigb.nc")
    sig_diwasp("11261sigb.nc")


def hobo_raw(glob_att, config_yaml):
    run_script("runhobocsv2cdf.py", glob_att, config_yaml)


def hobo_nc(nc_file):
    run_script("runhobocdf2nc.py", nc_file)


def test_hobo():
    hobo_raw("glob_att1168_hobowl.txt", "1168hwl_config.yaml")
    hobo_nc("11681hwl-raw.cdf")
    hobo_raw("glob_att1171_hobowl_baro.txt", "1171hwl_baro_config.yaml")
    hobo_nc("11711hwlb-raw.cdf")


def lisst_raw(glob_att, config_yaml):
    run_script("runlisstcsv2cdf.py", glob_att, config_yaml)


def lisst_nc(nc_file):
    run_script("runlisstcdf2nc.py", nc_file)


def test_lisst():
    lisst_raw("gatts_lisst_L0221705.txt", "config_lisst_L0221705.yaml")
    lisst_nc("lisst_L0221705-raw.cdf")


def tcm_raw(glob_att, config_yaml):
    run_script("runtcmcsv2cdf.py", glob_att, config_yaml)


def tcm_nc(nc_file):
    run_script("runtcmcdf2nc.py", nc_file)


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


def mc_raw(glob_att, config_yaml):
    run_script("runmcasc2cdf.py", glob_att, config_yaml)


def mc_nc(nc_file):
    run_script("runmccdf2nc.py", nc_file)


def test_mc():
    mc_raw("glob_att1126_mc.txt", "11263mc_config.yaml")
    mc_nc("11263mc-raw.cdf")


def sg_raw(glob_att, config_yaml):
    run_script("runots.py", "sgtid", "tid2cdf", glob_att, config_yaml)


def sg_nc(nc_file, atmpres, salwtemp):
    run_script(
        "runots.py",
        "sgtid",
        "cdf2nc",
        nc_file,
        "--atmpres",
        atmpres,
        "--salwtemp",
        salwtemp,
    )


def test_sg():
    sg_raw("sg_glob_att1126.txt", "11264sg_config.yaml")
    sg_nc("11264sg-tide-raw.cdf", "11264sg-atmpres.cdf", "11263mc-a.nc")


def sg_wv_raw(glob_att, config_yaml):
    run_script("runots.py", "sgwvs", "wb2cdf", glob_att, config_yaml)


def sg_wv_nc(nc_file, atmpres):
    run_script("runots.py", "sgwvs", "cdf2nc", nc_file, "--atmpres", atmpres)


def sg_wv_wvs(nc_file):
    run_script("runots.py", "sgwvs", "nc2waves", nc_file)


def test_sg_wvs():
    sg_wv_raw("sg_glob_att1126.txt", "11264sg_config.yaml")
    sg_wv_nc("11264sg-waves-raw.cdf", "11264sg-atmpres.cdf")
    sg_wv_wvs("11264sgb.nc")


def tb_raw(glob_att, config_yaml):
    run_script("runots.py", "tb", "csv2cdf", glob_att, config_yaml)


def tb_nc(nc_file, atmpres=None):
    if atmpres is not None:
        run_script("runots.py", "tb", "cdf2nc", nc_file, "--atmpres", atmpres)
    else:
        run_script("runots.py", "tb", "cdf2nc", nc_file)


def tb_wvs(nc_file):
    run_script("runots.py", "tb", "nc2waves", nc_file)


def test_tb():
    with zipfile.ZipFile("stglib/tests/data/example_TruBlue.zip", "r") as zip_ref:
        zip_ref.extractall("stglib/tests/data/")
    tb_raw("TB_glob_att.txt", "TB_config.yaml")
    tb_nc("example_TruBlue-raw.cdf")
    tb_wvs("example_TruBlue-cont-cal.nc")


def glx_dat(glob_att, config_yaml):
    run_script("runots.py", "glx", "dat2cdf", glob_att, config_yaml)


def glx_nc(nc_file, atmpres=None):
    run_script("runots.py", "glx", "cdf2nc", nc_file)


def glx_wvs(nc_file):
    run_script("runots.py", "glx", "nc2waves", nc_file)


def test_glx():
    glx_dat("glob_att_geolux_x600_202403.txt", "geolux_x600_202403_config.yaml")
    glx_nc("FRFx600_202403glx-raw.cdf")
    glx_wvs("FRFx600_202403glxb.nc")


def son_raw(glob_att, config_yaml):
    run_script("runots.py", "son", "raw2cdf", glob_att, config_yaml)


def son_nc(nc_file, height=None):
    if height is not None:
        run_script("runots.py", "son", "cdf2nc", nc_file, "--height", height)
    else:
        run_script("runots.py", "son", "cdf2nc", nc_file)


def son_xy(nc_file):
    run_script("runots.py", "son", "nc2xy", nc_file)


def test_son():
    son_raw("glob_att1126son.txt", "11265son_config.yaml")
    son_nc("11265son_5m-raw.cdf", "1126abs910s-cal.nc")
    son_xy("11265sonb_5m-a.nc")
