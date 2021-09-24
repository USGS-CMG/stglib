import numpy as np
import pandas as pd
import xarray as xr
import warnings

from .core import utils


def csv_to_cdf(metadata):
    """
    Process Aqua TROLL .csv file to a raw .cdf file
    """

    basefile = metadata["basefile"]

    try:
        df = read_aquatroll(
            basefile + ".csv",
            skiprows=metadata["skiprows"],
            encoding="utf-8",
            skipfooter=metadata["skipfooter"],
        )
    except UnicodeDecodeError:
        # try reading as Mac OS Western for old versions of Mac Excel
        df = read_aquatroll(
            basefile + ".csv",
            skiprows=metadata["skiprows"],
            encoding="mac-roman",
            skipfooter=metadata["skipfooter"],
        )

    ds = df_to_ds(df)

    metadata.pop("skiprows")
    metadata.pop("skipfooter")

    # write out metadata first, then deal exclusively with xarray attrs
    ds = utils.write_metadata(ds, metadata)

    del metadata

    ds = troll_shift_time(ds)

    if not utils.is_cf(ds):
        ds = utils.create_epic_times(ds)

    # configure file
    cdf_filename = ds.attrs["filename"] + "-raw.cdf"

    ds.to_netcdf(cdf_filename, unlimited_dims=["time"])

    print("Finished writing data to %s" % cdf_filename)

    return ds


def cdf_to_nc(cdf_filename):
    """
    Load a "raw" .cdf file and generate a processed .nc file
    """

    # Load raw .cdf data
    ds = xr.load_dataset(cdf_filename)

    # Clip data to in/out water times or via good_ens
    ds = utils.clip_ds(ds)

    # ds = ds_rename_vars(ds)
    #
    #
    # if "drop_vars" in ds.attrs:
    #     for k in ds.attrs["drop_vars"]:
    #         if k in ds:
    #             ds = ds.drop(k)
    #

    ds = compute_depth(ds)

    if "elev_offset" in ds.attrs:
        ds["waterlevel"] = ds["depth"] + ds.attrs["elev_offset"]

    ds = ds_add_attrs(ds)

    # ds = exo_qaqc(ds)
    #
    # assign min/max:
    ds = utils.add_min_max(ds)

    ds = utils.add_start_stop_time(ds)

    ds = add_delta_t(ds)

    # add lat/lon coordinates
    ds = utils.ds_add_lat_lon(ds)
    #
    # ds = ds_add_attrs(ds)
    #
    # # ds = utils.create_water_depth(ds)
    # ds = utils.create_nominal_instrument_depth(ds)
    #
    # ds = utils.no_p_create_depth(ds)
    #
    # # add lat/lon coordinates to each variable
    # for var in ds.variables:
    #     if (var not in ds.coords) and ("time" not in var):
    #         ds = utils.add_lat_lon(ds, var)
    #         ds = utils.no_p_add_depth(ds, var)
    #         # cast as float32
    #         ds = utils.set_var_dtype(ds, var)

    # Write to .nc file
    print("Writing cleaned/trimmed data to .nc file")
    nc_filename = ds.attrs["filename"] + "-a.nc"

    ds.to_netcdf(nc_filename, unlimited_dims=["time"])
    print("Done writing netCDF file", nc_filename)


def read_aquatroll(filnam, skiprows=69, encoding="utf-8", skipfooter=0):

    # so we don't get fallback parser warnings when specifying skipfooter
    warnings.filterwarnings("ignore", category=pd.errors.ParserWarning)

    df = pd.read_csv(
        filnam,
        skiprows=skiprows,
        skipfooter=skipfooter,
        infer_datetime_format=True,
        parse_dates=[0],
        encoding=encoding,
    )

    df.columns = df.columns.str.strip()

    md = get_metadata(filnam, encoding=encoding)

    df.attrs["type"] = md["ss"]
    df.attrs["sample_interval"] = f"{md['si']} {md['siu']}"
    df.attrs["samples_averaged"] = md["sa"]
    df.attrs["serial_number"] = md["sn"]
    df.attrs["INST_TYPE"] = md["de"]

    df.rename(
        columns={
            "Temperature (C)": "temperature",
            "Actual Conductivity (µS/cm)": "conductivity",
            "Pressure (kPa)": "pressure",
            "Pressure (PSI)": "pressure",
            "Date and Time (UTC)": "time",
            "Date and Time": "time",
        },
        inplace=True,
    )

    return df


def read_aquatroll_header(filnam, encoding="utf-8"):
    with open(filnam, encoding=encoding) as f:
        for line in f.readlines():
            if "Time Zone:" in line:
                # remove commas and only return the value,
                # not the 'Time Zone: ' part
                return line.replace(",", "").strip()[11:]


def compute_depth(ds):
    ds["salinity"] = compute_S(ds["temperature"], ds["conductivity"])
    ds["density"] = compute_density(ds["temperature"], ds["salinity"])
    ds["g"] = xr.ones_like(ds["density"]) * compute_g(np.deg2rad(48), 0)
    ds["depth"] = ds["pressure"] / ds["density"] / ds["g"]
    return ds


def troll_shift_time(ds):
    # remove time offsets which appear to be caused by errors in the sensor
    for second in [1, 2, 5, 9, 15, 45]:
        change = ds["time"].dt.second == second
        ds["time"].values[change] = ds["time"].values[change] - pd.Timedelta(
            f"{second}s"
        )

    # once we have this, see if this is a linear average sampling and if so, shift by that value divided by two
    if ds.attrs["type"] == "Linear Average":
        if ds.attrs["sample_interval"].split(" ")[1] != "secs":
            raise NotImplementedError(
                f"Can only shift time by seconds, not {ds.attrs['sample_interval'].split(' ')[1]}"
            )
        toshift = (
            ds.attrs["samples_averaged"]
            / float(ds.attrs["sample_interval"].split(" ")[0])
            / 2
        )
        ds = utils.shift_time(ds, toshift)

    return ds


def get_metadata(filnam, encoding="utf-8"):
    md = {}
    md["sn"] = 0
    md["ss"] = ""
    md["si"] = 0
    md["siu"] = ""
    md["sa"] = 0
    with open(filnam, encoding=encoding) as f:
        for line in f.readlines():
            if "Device," in line:
                cleanline = line.rstrip().split(",")
                md["de"] = cleanline[1]
            if "Serial Number," in line:
                cleanline = line.rstrip().split(",")
                md["sn"] = cleanline[1]
            if "Type," in line:
                cleanline = line.rstrip().split(",")
                md["ss"] = cleanline[2]
            if "Sample Interval," in line:
                cleanline = line.rstrip().split(",")
                md["si"] = float(cleanline[2])
                md["siu"] = cleanline[3]
            if "Samples Averaged," in line:
                cleanline = line.rstrip().split(",")
                md["sa"] = float(cleanline[2])

    return md


def df_to_ds(df):
    df = df.set_index("time")
    print(df.attrs)

    ds = df.to_xarray()
    for k in df.attrs:
        ds.attrs[k] = df.attrs[k]

    for k in ds:
        if "Unnamed" in k:
            ds = ds.drop(k)
        if "Seconds" in k:
            ds = ds.drop(k)

    return ds


def add_delta_t(ds):
    deltat = np.asscalar((ds["time"][1] - ds["time"][0]) / np.timedelta64(1, "s"))
    if not deltat.is_integer():
        warnings.warn("DELTA_T is not an integer; casting as int in attrs")

    ds.attrs["DELTA_T"] = int(deltat)

    return ds


def ds_add_attrs(ds):
    ds.attrs["institution"] = "U.S. Geological Survey"

    ds["time"].attrs["standard_name"] = "time"
    ds["time"].attrs["axis"] = "T"
    ds["time"].encoding = {"dtype": "int32"}

    ds["temperature"].attrs["units"] = "degree_C"
    ds["temperature"].attrs["standard_name"] = "sea_water_temperature"

    ds["conductivity"].attrs["standard_name"] = "sea_water_electrical_conductivity"
    ds["conductivity"].attrs["units"] = "uS cm-1"

    ds["salinity"].attrs["standard_name"] = "sea_water_practical_salinity"
    ds["salinity"].attrs["units"] = "1"

    ds["pressure"].attrs["long_name"] = "Pressure of water at pressure transducer"
    ds["pressure"].attrs["units"] = "kilopascal"

    ds["waterlevel"].attrs["long_name"] = "Water level relative to NAVD88"
    ds["waterlevel"].attrs["units"] = "m"
    ds["waterlevel"].attrs[
        "standard_name"
    ] = "sea_surface_height_above_geopotential_datum"

    if "g" in ds:
        ds["g"].attrs["long_name"] = "Gravitational acceleration"
        ds["g"].attrs["units"] = "m s-2"

    if "density" in ds:
        ds["density"].attrs["standard_name"] = "sea_water_density"
        ds["density"].attrs["units"] = "kg m-3"

    if "depth" in ds:
        ds["depth"].attrs["standard_name"] = "depth"
        ds["depth"].attrs[
            "long_name"
        ] = "Depth of pressure transducer below water surface"
        ds["depth"].attrs["units"] = "m"
        ds["depth"].attrs["positive"] = "down"

    return ds


def compute_g(φ, H):
    return (
        9.780356 * (1 + 0.0052885 * np.sin(φ) ** 2 - 0.0000059 * np.sin(2 * φ) ** 2)
        - 0.003086 * H
    )


def compute_Rt(T, AC):
    r0 = 29752.63
    r2 = 3.429338
    r1 = 830.5102
    r3 = -0.02193934
    return AC / (r0 + r1 * T + r2 * T ** 2 + r3 * T ** 3)


def compute_X(Rt):
    return 400 * Rt


def compute_Y(Rt):
    return 100 * Rt


def f(T):
    return (T - 15) / (1 + 0.0162 * (T - 15))


def compute_S(T, AC):
    a0 = 0.0080
    b0 = 0.0005
    a1 = -0.1692
    b1 = -0.0056
    a2 = 25.3851
    b2 = -0.0066
    a3 = 14.0941
    b3 = -0.0375
    a4 = -7.0261
    b4 = 0.0636
    a5 = 2.7081
    b5 = -0.0144

    Rt = compute_Rt(T, AC)
    X = compute_X(Rt)
    Y = compute_Y(Rt)

    return (
        a0
        + a1 * Rt ** (1 / 2)
        + a2 * Rt
        + a3 * Rt ** (3 / 2)
        + a4 * Rt ** 2
        + a5 * Rt ** (5 / 2)
    )
    +f(T) * (
        b0
        + b1 * Rt ** (1 / 2)
        + b2 * Rt
        + b3 * Rt ** (3 / 2)
        + b4 * Rt ** 2
        + b5 * Rt ** (5 / 2)
    )
    -a0 / (1 + 1.5 * X + X ** 2)
    -b0 * f(T) / (1 + Y ** (1 / 2) + Y ** (3 / 2))


def compute_density(T, S):
    ρ0 = (
        999.842594
        + 0.06793952 * T
        - 0.00909529 * T ** 2
        + 1.001685e-4 * T ** 3
        - 1.120083e-6 * T ** 4
        + 6.536332e-9 * T ** 5
    )
    a = (
        0.824493
        - 0.004089 * T
        + 7.6438e-5 * T ** 2
        - 8.2467e-7 * T ** 3
        + 5.3875e-9 * T ** 4
    )
    b = -0.00572466 + 1.0227e-4 * T - 1.6546e-6 * T ** 2
    c = 0.000483140

    return (ρ0 + a * S + b * S ** (3 / 2) + c * S ** 2) / 1000
