import numpy as np
import xarray as xr
from tqdm import tqdm

from ..core import utils, waves


def nc_to_waves(nc_filename):
    """
    Process burst data to wave statistics
    """

    ds = xr.load_dataset(nc_filename)

    spec = waves.make_waves_ds(ds)

    for k in ["wp_peak", "wh_4061", "wp_4060", "pspec"]:
        ds[k] = spec[k]

    dopuv = False
    if "puv" in ds.attrs:
        if ds.attrs["puv"].lower() == "true":
            dopuv = True

    if dopuv:
        ds = ds_puv(ds)

    # keep burst mean P_1 and P_1ac for reference
    for k in ["P_1", "P_1ac"]:
        if k in ds:
            ds[k] = ds[k].mean(dim="sample", keep_attrs=True)

    for k in [
        "burst",
        "sample",
        "u_1205",
        "v_1206",
        "w_1204",
        "vel1_1277",
        "vel2_1278",
        "vel3_1279",
        "AGC1_1221",
        "AGC2_1222",
        "AGC3_1223",
        "SNR1",
        "SNR2",
        "SNR3",
        "cor1_1285",
        "cor2_1286",
        "cor3_1287",
        "AnalogInput1",
        "AnalogInput2",
        "Hdg_1215",
        "Ptch_1216",
        "Roll_1217",
        "Bat_106",
        "Tx_1211",
        "SV_80",
    ]:
        if k in ds:
            ds = ds.drop(k)

    ds = utils.trim_max_wp(ds)

    ds = utils.trim_min_wh(ds)

    ds = utils.trim_max_wh(ds)

    ds = utils.trim_wp_ratio(ds)

    if dopuv:
        ds = waves.puv_qaqc(ds)

    # Add attrs
    ds = utils.ds_add_wave_attrs(ds)

    # assign min/max (need to do this after trimming):
    ds = utils.add_min_max(ds)

    nc_filename = ds.attrs["filename"] + "s-a.nc"

    ds.to_netcdf(nc_filename, unlimited_dims=["time"])
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])
    print("Done writing netCDF file", nc_filename)

    return ds


def ds_puv(ds):
    print("Running puv_quick")

    N, M = np.shape(ds["u_1205"].squeeze())

    desc = {
        "Hrmsp": "Hrms (=Hmo) from pressure",
        "Hrmsu": "Hrms from u,v",
        "ubr": "Representative orbital velocity amplitude in freq. band ( first_frequency_cutoff <= f <= last_frequency_cutoff ) (m/s)",
        "omegar": "Representative orbital velocity (radian frequency)",
        "Tr": "Representative orbital velocity period (s)",
        "Tpp": "Peak period from pressure (s)",
        "Tpu": "Peak period from velocity (s)",
        "phir": "Representative orbital velocity direction (angles from x-axis, positive ccw)",
        "azr": "Representative orb. velocity direction (deg; geographic azimuth; ambiguous =/- 180 degrees)",
        "ublo": "ubr in freq. band (f <= first_frequency_cutoff) (m/s)",
        "ubhi": "ubr in freq. band (f >= last_frequency_cutoff) (m/s)",
        "ubig": "ubr in infra-gravity freq. band (first_frequency_cutoff f <= 1/20) (m/s)",
        "Hrmsp_tail": "Hrms (=Hmo) from pressure with tail applied",
        "Hrmsu_tail": "Hrms from u,v with tail applied",
        "phir_tail": "Representative orbital velocity direction (angles from x-axis, positive ccw) with tail applied",
        "azr_tail": "Representative orb. velocity direction (deg; geographic azimuth; ambiguous =/- 180 degrees) with tail applied",
    }
    standard_name = {
        # "Hrmsp": "sea_surface_wave_significant_height",
        # "Hrmsu": "sea_surface_wave_significant_height",
        # "ubr": "Representative orbital velocity amplitude in freq. band ( first_frequency_cutoff <= f <= last_frequency_cutoff ) (m/s)",
        # "omegar": "Representative orbital velocity (radian frequency)",
        # "Tr": "Representative orbital velocity period (s)",
        "Tpp": "sea_surface_wave_period_at_variance_spectral_density_maximum",
        "Tpu": "sea_surface_wave_period_at_variance_spectral_density_maximum",
        "phir": "sea_surface_wave_from_direction",
        "azr": "sea_surface_wave_from_direction",
        # "ublo": "ubr in freq. band (f <= first_frequency_cutoff) (m/s)",
        # "ubhi": "ubr in freq. band (f >= last_frequency_cutoff) (m/s)",
        # "ubig": "ubr in infra-gravity freq. band (first_frequency_cutoff f <= 1/20) (m/s)",
        # "Hrmsp_tail": "sea_surface_wave_significant_height",
        # "Hrmsu_tail": "sea_surface_wave_significant_height",
        "phir_tail": "sea_surface_wave_from_direction",
        "azr_tail": "sea_surface_wave_from_direction",
    }

    puvs = {k: np.full_like(ds["time"].values, np.nan, dtype=float) for k in desc}

    if "P_1ac" in ds:
        pvar = "P_1ac"
    else:
        pvar = "P_1"

    if "puv_first_frequency_cutoff" in ds.attrs:
        first_frequency_cutoff = ds.attrs["puv_first_frequency_cutoff"]
    else:
        first_frequency_cutoff = 1 / 10

    if "puv_last_frequency_cutoff" in ds.attrs:
        last_frequency_cutoff = ds.attrs["puv_last_frequency_cutoff"]
    else:
        last_frequency_cutoff = 1 / 2.5

    for n in tqdm(range(N)):
        puv = waves.puv_quick(
            ds[pvar][n, :].values,
            ds["u_1205"][n, :],
            ds["v_1206"][n, :],
            ds[pvar][n, :].mean().values + ds.attrs["pressure_sensor_height"],
            ds.attrs["pressure_sensor_height"],
            ds.attrs["velocity_sample_volume_height"],
            1 / ds.attrs["sample_interval"],
            first_frequency_cutoff=first_frequency_cutoff,
            last_frequency_cutoff=last_frequency_cutoff,
        )

        for k in puvs:
            puvs[k][n] = puv[k]

    for k in puvs:
        ds["puv_" + k] = xr.DataArray(puvs[k], dims="time")
        if k in desc:
            ds["puv_" + k].attrs["description"] = desc[k]
        if k in standard_name:
            ds["puv_" + k].attrs["standard_name"] = standard_name[k]

    # add in Hs
    ds["puv_Hsp"] = np.sqrt(2) * ds["puv_Hrmsp"]
    ds["puv_Hsp"].attrs["description"] = "Hs computed via sqrt(2) * Hrmsp"
    ds["puv_Hsp"].attrs["standard_name"] = "sea_surface_wave_significant_height"
    ds["puv_Hsp"].attrs["units"] = "m"

    ds["puv_Hsu"] = np.sqrt(2) * ds["puv_Hrmsu"]
    ds["puv_Hsu"].attrs["description"] = "Hs computed via sqrt(2) * Hrmsu"
    ds["puv_Hsu"].attrs["standard_name"] = "sea_surface_wave_significant_height"
    ds["puv_Hsu"].attrs["units"] = "m"

    ds["puv_Hsp_tail"] = np.sqrt(2) * ds["puv_Hrmsp_tail"]
    ds["puv_Hsp_tail"].attrs[
        "description"
    ] = "Hs computed via sqrt(2) * Hrmsp with tail applied"
    ds["puv_Hsp_tail"].attrs["standard_name"] = "sea_surface_wave_significant_height"
    ds["puv_Hsp_tail"].attrs["units"] = "m"

    ds["puv_Hsu_tail"] = np.sqrt(2) * ds["puv_Hrmsu_tail"]
    ds["puv_Hsu_tail"].attrs[
        "description"
    ] = "Hs computed via sqrt(2) * Hrmsu with tail applied"
    ds["puv_Hsu_tail"].attrs["standard_name"] = "sea_surface_wave_significant_height"
    ds["puv_Hsu_tail"].attrs["units"] = "m"

    ds["puv_Tpp"].attrs["units"] = "s"

    ds["puv_Tpu"].attrs["units"] = "s"

    ds["puv_azr"].attrs["units"] = "degrees"
    ds["puv_azr_tail"].attrs["units"] = "degrees"

    return ds
