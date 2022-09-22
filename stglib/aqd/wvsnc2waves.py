import xarray as xr
import numpy as np

from ..core import utils, waves
from . import aqdutils
from tqdm import tqdm


def nc_to_waves(nc_filename):

    ds = xr.load_dataset(nc_filename, decode_times=False)

    if utils.is_cf(ds):
        for k in ds:
            if "_time" in k:
                ds = ds.drop(k)
        ds = xr.decode_cf(ds)
    else:
        ds = utils.epic_to_cf_time(ds)
        ds = utils.create_epic_times(ds)

    spec = waves.make_waves_ds(ds)

    for k in ["wp_peak", "wh_4061", "wp_4060", "pspec"]:
        ds[k] = spec[k]

    dopuv = True
    if dopuv:
        print("Computing PUV")
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
        if ds.attrs["orientation"] != "UP":
            raise NotImplementedError(
                "PUV currently only works on UP oriented instruments"
            )
        # only works on UP because we assume T = T_orig
        # need to pass .values because passing the DataArray is MUCH slower
        u, v, w = aqdutils.coord_transform(
            ds["vel1_1277"].squeeze().values / 1000,
            ds["vel2_1278"].squeeze().values / 1000,
            ds["vel3_1279"].squeeze().values / 1000,
            ds["Hdg_1215"].squeeze().values,
            ds["Ptch_1216"].squeeze().values,
            ds["Roll_1217"].squeeze().values,
            ds["TransMatrix"].squeeze().values,
            ds["TransMatrix"].squeeze().values,
            ds.attrs["AQDCoordinateSystem"],
        )
        ds["u_1205"] = xr.DataArray(u, dims=("time", "sample"))
        ds["v_1206"] = xr.DataArray(v, dims=("time", "sample"))
        ds["w_1204"] = xr.DataArray(w, dims=("time", "sample"))
        N, M = np.shape(ds["vel1_1277"].squeeze())
        puvs = {
            k: np.full_like(ds["time"].values, np.nan, dtype=float)
            for k in desc
            # "Hrmsp": np.full_like(ds["time"].values, np.nan, dtype=float),
            # "Hrmsu": np.full_like(ds["time"].values, np.nan, dtype=float),
            # "ubr": np.full_like(ds["time"].values, np.nan, dtype=float),
            # "omegar": np.full_like(ds["time"].values, np.nan, dtype=float),
            # "Tr": np.full_like(ds["time"].values, np.nan, dtype=float),
            # "Tpp": np.full_like(ds["time"].values, np.nan, dtype=float),
            # "Tpu": np.full_like(ds["time"].values, np.nan, dtype=float),
            # "phir": np.full_like(ds["time"].values, np.nan, dtype=float),
            # "azr": np.full_like(ds["time"].values, np.nan, dtype=float),
            # "ublo": np.full_like(ds["time"].values, np.nan, dtype=float),
            # "ubhi": np.full_like(ds["time"].values, np.nan, dtype=float),
            # "ubig": np.full_like(ds["time"].values, np.nan, dtype=float),
            # "Hrmsp_tail": np.full_like(ds["time"].values, np.nan, dtype=float),
            # "Hrmsu_tail": np.full_like(ds["time"].values, np.nan, dtype=float),
        }

        print("Running puv_quick")
        for n in tqdm(range(N)):
            try:
                puv = waves.puv_quick(
                    ds["P_1ac"][n, :].values,
                    u[n, :],
                    v[n, :],
                    ds["P_1ac"][n, :].mean().values
                    + ds.attrs["initial_instrument_height"],
                    ds.attrs["initial_instrument_height"],
                    ds.attrs["initial_instrument_height"]
                    + ds.attrs["center_first_bin"],
                    1 / ds.attrs["sample_interval"],
                    first_frequency_cutoff=1 / 10,
                    last_frequency_cutoff=1 / 2.5,
                )
                for k in puvs:
                    puvs[k][n] = puv[k]

            except AttributeError:  # puv_quick will fail on some values and return AttributeError: 'float' object has no attribute 'astype'
                continue

        for k in puvs:
            ds["puv_" + k] = xr.DataArray(puvs[k], dims="time")
            ds["puv_" + k].attrs["description"] = desc[k]
            ds["puv_" + k].attrs["standard_name"] = standard_name[k]

        # add in Hs
        ds["puv_Hsp"] = np.sqrt(2) * ds["puv_Hrmsp"]
        ds["puv_Hsp"].attrs["description"] = "Hs computed via sqrt(2) * Hrmsp"
        ds["puv_Hsp"].attrs["standard_name"] = "sea_surface_wave_significant_height"
        ds["puv_Hsu"] = np.sqrt(2) * ds["puv_Hrmsu"]
        ds["puv_Hsu"].attrs["description"] = "Hs computed via sqrt(2) * Hrmsu"
        ds["puv_Hsu"].attrs["standard_name"] = "sea_surface_wave_significant_height"
        ds["puv_Hsp_tail"] = np.sqrt(2) * ds["puv_Hrmsp_tail"]
        ds["puv_Hsp_tail"].attrs[
            "description"
        ] = "Hs computed via sqrt(2) * Hrmsp with tail applied"
        ds["puv_Hsp_tail"].attrs[
            "standard_name"
        ] = "sea_surface_wave_significant_height"
        ds["puv_Hsu_tail"] = np.sqrt(2) * ds["puv_Hrmsu_tail"]
        ds["puv_Hsu_tail"].attrs[
            "description"
        ] = "Hs computed via sqrt(2) * Hrmsu with tail applied"
        ds["puv_Hsu_tail"].attrs[
            "standard_name"
        ] = "sea_surface_wave_significant_height"

    # ds = utils.create_water_depth(ds)

    # Remove old variables as we just want to keep the wave statistics
    keys = [
        "P_1",
        "P_1ac",
        "sample",
        "Tx_1211",
        "vel1_1277",
        "vel2_1278",
        "vel3_1279",
        "U",
        "V",
        "W",
        "avgamp1",
        "avgamp2",
        "avgamp3",
        "AGC1_1221",
        "AGC2_1222",
        "AGC3_1223",
        "TransMatrix",
        "nrecs",
        "burst",
        "soundspeed",
        "Battery",
        "Hdg_1215",
        "Ptch_1216",
        "Roll_1217",
    ]

    if dopuv:
        keys.remove("sample")

    for k in keys:
        if k in ds:
            ds = ds.drop(k)

    ds = utils.trim_max_wp(ds)

    ds = utils.trim_min_wh(ds)

    ds = utils.trim_max_wh(ds)

    ds = utils.trim_wp_ratio(ds)

    ds = waves.puv_qaqc(ds)

    # Add attrs
    ds = utils.ds_add_attrs(ds)

    nc_filename = ds.attrs["filename"] + "wvs-a.nc"

    ds.to_netcdf(nc_filename, unlimited_dims=["time"])
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])

    print("Done creating", nc_filename)

    return ds
