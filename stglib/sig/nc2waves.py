import re

import numpy as np
import xarray as xr

from ..core import qaqc, utils, waves


def nc_to_waves(nc_filename):
    """
    Process burst data to wave statistics
    """
    ds = xr.open_dataset(nc_filename)

    # check to see if need to make wave burst from continuous data
    if (ds.attrs["sample_mode"] == "CONTINUOUS") and ("wave_interval" in ds.attrs):
        # check for wave_start_time attrs
        if "wave_start_time" in ds.attrs:
            print(
                f"trimming continuous data to start at users specified wave start time {ds.attrs['wave_start_time']}"
            )
            ds = ds.sel(
                time=slice(np.datetime64(ds.attrs["wave_start_time"]), ds["time"][-1])
            )
        # make wave burst ncfile from continuous data if wave_interval is specified
        ds = make_wave_bursts(ds)

    spec = waves.make_waves_ds(ds)
    ds = utils.ds_add_waves_history(ds)

    for k in ["wp_peak", "wh_4061", "wp_4060", "pspec"]:
        ds[k] = spec[k]

    dopuv = False
    if "puv" in ds.attrs:
        if ds.attrs["puv"].lower() == "true":
            dopuv = True

    if dopuv:
        ds = do_puv(ds)

    ds = utils.create_water_depth_var(ds)

    for k in [
        "burst",
        "sample",
        "P_1",
        "P_1ac",
        "u_1205",
        "v_1206",
        "w_1204",
        "Tx_1211",
        "SV_80",
        "vel",
        "amp",
        "cor",
        "brangeAST",
        "TransMatrix",
    ]:
        if k in ds:
            ds = ds.drop_vars(k)

    # ds = qaqc.drop_vars(ds)

    ds = drop_unused_dims(ds)

    ds = utils.trim_max_wp(ds)

    ds = utils.trim_min_wh(ds)

    ds = utils.trim_max_wh(ds)

    ds = utils.trim_wp_ratio(ds)

    if dopuv:
        ds = waves.puv_qaqc(ds)

    # Add attrs
    ds = utils.ds_add_wave_attrs(ds)

    # clean-up attrs
    ds = drop_attrs(ds)

    # assign min/max (need to do this after trimming):
    ds = utils.add_min_max(ds)

    # round time to minutes to make sure fits in dtype i4. will be fine for wave burst start times
    ds["time"] = ds["time"].dt.round("min")

    if dopuv:
        nc_filename = ds.attrs["filename"] + "s_puvq-a.nc"
    else:
        nc_filename = ds.attrs["filename"] + "s-a.nc"

    ds.to_netcdf(nc_filename, unlimited_dims=["time"])
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])

    print("Done writing netCDF file", nc_filename)

    return ds


def nc_to_diwasp(nc_filename):
    """
    Process burst data to make DIWASP derived wave statistics and spectra
    """

    ds = xr.open_dataset(nc_filename)

    # check to see if need to make wave burst from continuous data
    if (ds.attrs["sample_mode"] == "CONTINUOUS") and ("wave_interval" in ds.attrs):
        # check for wave_start_time attrs
        if "wave_start_time" in ds.attrs:
            print(
                f"trimming continuous data to start at users specified wave start time {ds.attrs['wave_start_time']}"
            )
            ds = ds.sel(
                time=slice(np.datetime64(ds.attrs["wave_start_time"]), ds["time"][-1])
            )
        # make wave burst ncfile from continuous data if wave_interval is specified
        ds = make_wave_bursts(ds)

    if "diwasp" not in ds.attrs:
        # if not specified choose 'suv' method for signature
        ds.attrs["diwasp"] = "suv"

    if "diwasp" in ds.attrs:
        data_type = ds.attrs["diwasp"]
        if data_type not in ["suv", "puv", "elev", "pres"]:
            raise ValueError(
                f"data type {data_type} is not recognized current options are ['suv', 'puv', 'elev', 'pres']"
            )

    if "diwasp_ibin" in ds.attrs:
        ibin = ds.attrs["diwasp_ibin"]
    else:
        ibin = 0

    if data_type == "puv" or data_type == "suv":
        print(f"Running DIWASP using {data_type} input data")
        layout = make_diwasp_layout(ds, data_type=data_type, ibin=ibin)
        diwasp = waves.make_diwasp_puv_suv(
            ds, data_type=data_type, layout=layout, ibin=ibin
        )
        ds = utils.ds_add_pydiwasp_history(ds)

    elif data_type == "elev" or data_type == "pres":
        print(f"Running DIWASP using {ds.attrs['diwasp']} input data")
        layout = make_diwasp_layout(ds, data_type=data_type)
        diwasp = waves.make_diwasp_elev_pres(ds, data_type=data_type, layout=layout)
        ds = utils.ds_add_pydiwasp_history(ds)
    else:
        raise ValueError(
            f"DIWASP input type {ds.attrs['diwasp']} is not currently supported for {ds.attrs['instument_type']} in stglib"
        )

    for k in [
        "diwasp_frequency",
        "diwasp_direction",
        "diwasp_hs",
        "diwasp_tp",
        "diwasp_tm",
        "diwasp_dtp",
        "diwasp_dp",
        "diwasp_dm",
        "diwasp_dspec",
        "diwasp_fspec",
    ]:

        if k in diwasp:
            ds[k] = diwasp[k]

    # add diwasp attrs
    for k in diwasp.attrs:
        ds.attrs[k] = diwasp.attrs[k]

    # rename Fspec based on input datatype
    ds = utils.rename_diwasp_fspec(ds)

    ds = utils.create_water_depth_var(ds)

    for k in [
        "burst",
        "sample",
        "P_1",
        "P_1ac",
        "u_1205",
        "v_1206",
        "w_1204",
        "Tx_1211",
        "SV_80",
        "vel",
        "amp",
        "cor",
        "brangeAST",
        "TransMatrix",
    ]:
        if k in ds:
            ds = ds.drop_vars(k)

    # ds = qaqc.drop_vars(ds)

    # use "epic" names
    ds.attrs["diwasp_names"] = "epic"
    ds = utils.rename_diwasp_wave_vars(ds)

    ds = drop_unused_dims(ds)

    ds = utils.trim_max_wp(ds)

    ds = utils.trim_min_wh(ds)

    ds = utils.trim_max_wh(ds)

    ds = utils.trim_wp_ratio(ds)

    # Add attrs
    ds = utils.ds_add_wave_attrs(ds)

    # clean-up attrs
    ds = drop_attrs(ds)

    # assign min/max (need to do this after trimming):
    ds = utils.add_min_max(ds)

    # round time to minutes to make sure fits in dtype i4. will be fine for wave burst start times
    ds["time"] = ds["time"].dt.round("min")

    nc_filename = ds.attrs["filename"] + "s_diwasp-a.nc"

    ds.to_netcdf(nc_filename, unlimited_dims=["time"])
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])

    print("Done writing netCDF file", nc_filename)

    return ds


def make_diwasp_layout(ds, data_type=None, ibin=None):
    """
    Make layout for DIWASP wave processing input
    """

    if data_type == "suv" or data_type == "puv":
        # datatypes = ["pres", "velx", "vely"]
        sxyz = [0, 0, ds.attrs["initial_instrument_height"]]
        if ds.attrs["orientation"].lower() == "up":
            velz = ds.attrs["initial_instrument_height"] + ds["bindist"][ibin].values

        elif ds.attrs["orientation"].lower() == "down":
            velz = ds.attrs["initial_instrument_height"] - ds["bindist"][ibin].values

        uxyz = vxyz = [0, 0, velz]

        layout = np.array([sxyz, uxyz, vxyz]).T

    elif data_type == "elev" or data_type == "pres":
        # datatypes=['elev']
        sxyz = np.atleast_2d([0, 0, ds.attrs["initial_instrument_height"]]).T
        layout = sxyz

    return layout


def make_wave_bursts(ds):
    # wave_interval is [sec] interval for wave statistics for continuous data
    ds.attrs["wave_samples_per_burst"] = int(
        ds.attrs["wave_interval"] / ds.attrs["sample_interval"]
    )

    r = np.shape(ds.P_1)[0]
    mod = r % ds.attrs["wave_samples_per_burst"]
    if mod:
        print(
            "Number of rows is not a multiple of samples_per_burst; truncating to last full burst"
        )
        ds = ds.sel(time=ds.time[0:-mod])

    ds.time.encoding.pop("dtype")

    ds["timenew"] = xr.DataArray(
        ds.time[0 :: int(ds.attrs["wave_samples_per_burst"])].values, dims="timenew"
    )

    ds["samplenew"] = xr.DataArray(
        range(ds.attrs["wave_samples_per_burst"]), dims="samplenew"
    )

    for v in ["P_1", "P_1ac", "Tx_1211", "brangeAST"]:
        if v in ds:
            attrsbak = ds[v].attrs
            ds[v] = xr.DataArray(
                np.reshape(ds[v].values, (-1, int(ds.attrs["wave_samples_per_burst"]))),
                dims=["timenew", "samplenew"],
            )
            ds[v].attrs = attrsbak

    for v in ["u_1205", "v_1206", "w_1204"]:
        if v in ds:
            attrsbak = ds[v].attrs
            # check dims
            if ds[v].dims == ("time", "z"):
                ds[v] = xr.DataArray(
                    np.reshape(
                        ds[v].values,
                        (-1, int(ds.attrs["wave_samples_per_burst"]), len(ds["z"])),
                    ).transpose(0, 2, 1),
                    dims=["timenew", "z", "samplenew"],
                )
                # u=np.reshape(ds['u_1205'].values, (-1,int(ds.attrs["samples_per_burst"]), len(ds['z']))).transpose(0,2,1)
                ds[v].attrs = attrsbak
            else:
                raise ValueError(
                    f"{v} dimensions {ds[v].dims} are not as required ('time','z') to shape into wave burst"
                )

    for v in ["vel", "cor", "amp"]:
        if v in ds:
            attrsbak = ds[v].attrs
            if ds[v].dims == ("beam", "time", "z"):
                ds[v] = xr.DataArray(
                    np.reshape(
                        ds[v].values,
                        (
                            len(ds["beam"]),
                            -1,
                            int(ds.attrs["wave_samples_per_burst"]),
                            len(ds["z"]),
                        ),
                    ).transpose(0, 1, 3, 2),
                    dims=["beam", "timenew", "z", "samplenew"],
                )
                # u=np.reshape(ds['u_1205'].values, (-1,int(ds.attrs["samples_per_burst"]), len(ds['z']))).transpose(0,2,1)
                ds[v].attrs = attrsbak
            else:
                raise ValueError(
                    f"{v} dimensions {ds[v].dims} are not as required ('beam','time','z') to shape into wave burst"
                )

    ds = ds.rename({"time": "timeold"})
    ds = ds.rename({"timenew": "time"})
    ds = ds.drop_dims("timeold")

    if "sample" in ds:
        ds = ds.rename({"sample": "sampleold"})
        ds = ds.rename({"samplenew": "sample"})
        ds = ds.drop_vars("sampleold")
    else:
        ds = ds.rename({"samplenew": "sample"})

    return ds


def drop_unused_dims(ds):
    """only keep dims that will be in the final files"""
    thedims = []
    for v in ds.data_vars:
        for x in ds[v].dims:
            thedims.append(x)

    for x in ds.dims:
        if x not in thedims:
            ds = ds.drop_vars(x)

    return ds


def drop_attrs(ds):
    """Drop some global attrs"""

    rm = []  # initialize list of attrs to be removed
    exclude = []  # initialize attrs to exclude

    for j in ds.attrs:
        if re.match("^SIG", j):
            rm.append(j)

    for k in rm:
        if k not in exclude:
            del ds.attrs[k]

    return ds


def do_puv(ds):
    print("Running puv_quick")

    for k in ["initial_instrument_height"]:
        if k not in ds.attrs:
            raise KeyError(f"{k} must be specified to run PUV")

    if "puv_bin" in ds.attrs:
        ibin = ds.attrs["puv_bin"]
    else:
        ibin = 0

    N, M = np.shape(ds["u_1205"].isel(z=ibin).squeeze())

    if "puv_first_frequency_cutoff" in ds.attrs:
        first_frequency_cutoff = ds.attrs["puv_first_frequency_cutoff"]
    else:
        first_frequency_cutoff = 1 / 33

    if "puv_last_frequency_cutoff" in ds.attrs:
        last_frequency_cutoff = ds.attrs["puv_last_frequency_cutoff"]
    else:
        last_frequency_cutoff = 1 / 3.3

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

    if ds.attrs["orientation"].lower() == "up":
        huv = ds.attrs["initial_instrument_height"] + ds["bindist"][ibin].values

    elif ds.attrs["orientation"].lower() == "down":
        huv = ds.attrs["initial_instrument_height"] - ds["bindist"][ibin].values

    puvs = waves.puv_quick_vectorized(
        ds[pvar].squeeze(),
        ds["u_1205"].isel(z=ibin).squeeze(),
        ds["v_1206"].isel(z=ibin).squeeze(),
        ds[pvar].squeeze().mean(dim="sample").values
        + ds.attrs["initial_instrument_height"],
        ds.attrs["initial_instrument_height"],
        huv,
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
