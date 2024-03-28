import numpy as np
import xarray as xr

from ..core import utils, waves


def nc_to_waves(nc_filename):
    """
    Process burst data to wave statistics
    """

    ds = xr.load_dataset(nc_filename)

    spec = waves.make_waves_ds(ds)

    for k in ["wp_peak", "wh_4061", "wp_4060", "pspec", "pspec_jonswap"]:
        ds[k] = spec[k]

    dopuv = False
    if "puv" in ds.attrs:
        if ds.attrs["puv"].lower() == "true":
            dopuv = True

    if dopuv:
        ds = do_puv(ds)

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


def do_puv(ds):
    print("Running puv_quick")

    for k in ["pressure_sensor_height", "velocity_sample_volume_height"]:
        if k not in ds.attrs:
            raise KeyError(f"{k} must be specified to run PUV")

    N, M = np.shape(ds["u_1205"].squeeze())

    if "puv_first_frequency_cutoff" in ds.attrs:
        first_frequency_cutoff = ds.attrs["puv_first_frequency_cutoff"]
    else:
        first_frequency_cutoff = 1 / 10

    if "puv_last_frequency_cutoff" in ds.attrs:
        last_frequency_cutoff = ds.attrs["puv_last_frequency_cutoff"]
    else:
        last_frequency_cutoff = 1 / 2.5

    desc = {
        "Hrmsp": f"Hrms (=Hmo) from pressure in freq. band {first_frequency_cutoff} <= f <= {last_frequency_cutoff}",
        "Hrmsu": f"Hrms from u,v in freq. band {first_frequency_cutoff} <= f <= {last_frequency_cutoff}",
        "ubr": f"Representative orbital velocity amplitude in freq. band {first_frequency_cutoff} <= f <= {last_frequency_cutoff}",
        "omegar": f"Representative orbital velocity (radian frequency) in freq. band {first_frequency_cutoff} <= f <= {last_frequency_cutoff}",
        "Tr": f"Representative orbital velocity period in freq. band {first_frequency_cutoff} <= f <= {last_frequency_cutoff}",
        "Tpp": f"Peak period from pressure in freq. band {first_frequency_cutoff} <= f <= {last_frequency_cutoff}",
        "Tpu": f"Peak period from velocity in freq. band {first_frequency_cutoff} <= f <= {last_frequency_cutoff}",
        "phir": "Representative orbital velocity direction (angles from x-axis, positive ccw)",
        "azr": "Representative orb. velocity direction (deg; geographic azimuth; ambiguous =/- 180 degrees)",
        "ublo": f"ubr in freq. band f <= {first_frequency_cutoff}",
        "ubhi": f"ubr in freq. band f >= {last_frequency_cutoff}",
        "ubig": f"ubr in infra-gravity freq. band {first_frequency_cutoff} f <= 1/20",
        "Hrmsp_tail": "Hrms (=Hmo) from pressure with f^-4 tail applied",
        "Hrmsu_tail": "Hrms from u,v with f^-4 tail applied",
        "phir_tail": "Representative orbital velocity direction (angles from x-axis, positive ccw) with f^-4 tail applied",
        "azr_tail": "Representative orb. velocity direction (deg; geographic azimuth; ambiguous =/- 180 degrees) with f^-4 tail applied",
        "Snp": f"Pressure-derived non-directional wave energy spectrum in freq. band {first_frequency_cutoff} <= f <= {last_frequency_cutoff}",
        "Snp_tail": "Pressure-derived non-directional wave energy spectrum with f^-4 tail applied",
        "Snu": f"Velocity-derived non-directional wave energy spectrum in freq. band {first_frequency_cutoff} <= f <= {last_frequency_cutoff}",
        "Snu_tail": "Velocity-derived non-directional wave energy spectrum with f^-4 tail applied",
        "frequencies": "Frequency",
        "fclip": "Frequency",
    }
    standard_name = {
        "Tpp": "sea_surface_wave_period_at_variance_spectral_density_maximum",
        "Tpu": "sea_surface_wave_period_at_variance_spectral_density_maximum",
        "phir": "sea_surface_wave_from_direction",
        "phir_tail": "sea_surface_wave_from_direction",
        "azr": "sea_surface_wave_from_direction",
        "azr_tail": "sea_surface_wave_from_direction",
        "Snp": "sea_surface_wave_variance_spectral_density",
        "Snp_tail": "sea_surface_wave_variance_spectral_density",
        "Snu": "sea_surface_wave_variance_spectral_density",
        "Snu_tail": "sea_surface_wave_variance_spectral_density",
    }
    unit = {
        "Hrmsp": "m",
        "Hrmsu": "m",
        "ubr": "m s-1",
        "omegar": "rad s-1",
        "Tr": "s",
        "Tpp": "s",
        "Tpu": "s",
        "phir": "radians",
        "phir_tail": "radians",
        "azr": "degrees",
        "azr_tail": "degrees",
        "ublo": "m s-1",
        "ubhi": "m s-1",
        "ubig": "m s-1",
        "Hrmsp_tail": "m",
        "Hrmsu_tail": "m",
        "Snp": "m2/Hz",
        "Snp_tail": "m2/Hz",
        "Snu": "m2/Hz",
        "Snu_tail": "m2/Hz",
    }

    # puvs = {k: np.full_like(ds["time"].values, np.nan, dtype=float) for k in desc}
    puvs = {k: [] for k in desc}

    if "P_1ac" in ds:
        pvar = "P_1ac"
    else:
        pvar = "P_1"

    puvs = waves.puv_quick_vectorized(
        ds[pvar].squeeze(),
        ds["u_1205"].squeeze(),
        ds["v_1206"].squeeze(),
        ds[pvar].squeeze().mean(dim="sample").values
        + ds.attrs["pressure_sensor_height"],
        ds.attrs["pressure_sensor_height"],
        ds.attrs["velocity_sample_volume_height"],
        1 / ds.attrs["sample_interval"],
        first_frequency_cutoff=first_frequency_cutoff,
        last_frequency_cutoff=last_frequency_cutoff,
    )

    ds["puv_frequency"] = xr.DataArray(
        puvs["frequencies"],
        dims="puv_frequency",
        attrs={"standard_name": "sea_surface_wave_frequency", "units": "Hz"},
    )
    ds["puv_frequency_clipped"] = xr.DataArray(
        puvs["fclip"],
        dims="puv_frequency_clipped",
        attrs={"standard_name": "sea_surface_wave_frequency", "units": "Hz"},
    )

    for k in puvs:
        if k == "frequencies" or k == "fclip":
            continue
        if puvs[k].ndim == 1:
            ds["puv_" + k] = xr.DataArray(puvs[k], dims="time")
        elif puvs[k].ndim == 2:
            try:
                ds["puv_" + k] = xr.DataArray(puvs[k], dims=["time", "puv_frequency"])
            except ValueError:
                ds["puv_" + k] = xr.DataArray(
                    puvs[k], dims=["time", "puv_frequency_clipped"]
                )
        if k in desc:
            ds["puv_" + k].attrs["description"] = desc[k]
        if k in standard_name:
            ds["puv_" + k].attrs["standard_name"] = standard_name[k]
        if k in unit:
            ds["puv_" + k].attrs["units"] = unit[k]

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
    ] = "Hs computed via sqrt(2) * Hrmsp with f^-4 tail applied"
    ds["puv_Hsp_tail"].attrs["standard_name"] = "sea_surface_wave_significant_height"
    ds["puv_Hsp_tail"].attrs["units"] = "m"

    ds["puv_Hsu_tail"] = np.sqrt(2) * ds["puv_Hrmsu_tail"]
    ds["puv_Hsu_tail"].attrs[
        "description"
    ] = "Hs computed via sqrt(2) * Hrmsu with f^-4 tail applied"
    ds["puv_Hsu_tail"].attrs["standard_name"] = "sea_surface_wave_significant_height"
    ds["puv_Hsu_tail"].attrs["units"] = "m"

    return ds
