import warnings

import numpy as np
import pandas as pd
import xarray as xr

from .core import qaqc, utils


def csv_to_cdf(metadata):
    """
    Process LISST .csv file to a raw .cdf file
    """

    basefile = metadata["basefile"]

    ds = read_lisst(basefile + ".csv")

    # write out metadata first, then deal exclusively with xarray attrs
    ds = utils.write_metadata(ds, metadata)

    del metadata

    ds = utils.ensure_cf(ds)

    ds = read_lop(ds)

    cdf_filename = ds.attrs["filename"] + "-raw.cdf"

    ds.to_netcdf(cdf_filename, unlimited_dims=["time"])

    print(f"Finished writing data to {cdf_filename}")

    return ds


def cdf_to_nc(cdf_filename, atmpres=False):
    """
    Load a raw .cdf file and generate a processed .nc file
    """

    # Load raw .cdf data
    ds = xr.load_dataset(cdf_filename)

    # remove units in case we change and we can use larger time steps
    ds.time.encoding.pop("units")

    # Infer if a burst deployment depending on whether there is a large DT
    ds = check_and_reshape_burst(ds)

    # Clip data to in/out water times or via good_ens
    ds = utils.clip_ds(ds)

    # QAQC
    ds = qaqc.call_qaqc(ds)

    ds = utils.add_min_max(ds, exclude_vars=["RSlower", "RSmedian", "RSupper"])

    ds = utils.add_start_stop_time(ds)

    ds = utils.add_delta_t(ds)

    ds = utils.ds_add_lat_lon(ds)

    # Write to .nc file
    print("Writing cleaned/trimmed data to .nc file")
    if "sample" in ds:
        nc_filename = ds.attrs["filename"] + "b-cal.nc"
    else:
        nc_filename = ds.attrs["filename"] + "-a.nc"

    ds.to_netcdf(
        nc_filename, unlimited_dims=["time"], encoding={"time": {"dtype": "i4"}}
    )
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])
    print("Done writing netCDF file", nc_filename)

    if "sample" in ds:
        print("Writing burst-averaged data to .nc file")
        nc_filename = ds.attrs["filename"] + "-a.nc"
        dsm = ds.mean(dim="sample", keep_attrs=True)
        # need to recompute min/max on averaged values
        dsm = utils.add_min_max(dsm, exclude_vars=["RSlower", "RSmedian", "RSupper"])

        dsm.to_netcdf(
            nc_filename, unlimited_dims=["time"], encoding={"time": {"dtype": "i4"}}
        )
        utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])
        print("Done writing netCDF file", nc_filename)


def read_lisst(f):
    """Read data from exported LISST CSV file"""
    vcs = [f"vc{n:02}" for n in range(1, 37)]
    cols = [
        "LaserTransmissionSensor",
        "SupplyVoltage",
        "AnalogInput1",
        "LaserReferenceSensor",
        "Depth",
        "Temperature",
        "Year",
        "Month",
        "Day",
        "Hour",
        "Minute",
        "Second",
        "AnalogInput2",
        "MeanDiameter",
        "TotalVolumeConcentration",
        "RelativeHumidity",
        "AccelerometerX",
        "AccelerometerY",
        "AccelerometerZ",
        "RawPressureMSB",
        "RawPressureLSBs",
        "AmbientLight",
        "AnalogInput3",
        "ComputedOpticalTransmissionOverPath",
        "BeamAttenuation",
    ]

    df = pd.read_csv(f, names=np.hstack([vcs, cols]), skipinitialspace=True)

    df["time"] = pd.to_datetime(
        df["Year"].astype(str)
        + "-"
        + df["Month"].astype(str)
        + "-"
        + df["Day"].astype(str)
        + " "
        + df["Hour"].astype(str)
        + ":"
        + df["Minute"].astype(str)
        + ":"
        + df["Second"].astype(str)
    )

    df = df.drop(
        columns=[
            "Year",
            "Month",
            "Day",
            "Hour",
            "Minute",
            "Second",
        ]
    )

    df = df.set_index("time")

    dfvc = df[[f"vc{n:02}" for n in range(1, 37)]]
    df = df.drop(columns=[f"vc{n:02}" for n in range(1, 37)])

    ds = df.to_xarray()
    ds["ring"] = xr.DataArray(np.arange(1, 37), dims="ring")
    ds["vc"] = xr.DataArray(dfvc.values, dims=("time", "ring"))

    ds["RSmedian"], ds["RSlower"], ds["RSupper"] = get_ringsizes()

    ds = set_metadata(ds)

    return ds


def set_metadata(ds):
    """
    Set metadata. Values obtained from:

    LISST-200X
    Particle Size Analyzer
    Including
    LISST-HAB & LISST-Black
    User’s Manual
    Version 2.3
    March, 2022

    Appendix C: Data File Formats (pp. 81)
    """

    ds["LaserTransmissionSensor"].attrs["units"] = "mW"
    ds["LaserTransmissionSensor"].attrs["long_name"] = "Laser transmission sensor"

    ds["SupplyVoltage"].attrs["units"] = "V"
    ds["SupplyVoltage"].attrs["long_name"] = "Supply voltage"

    ds["AnalogInput1"].attrs["units"] = "V"
    ds["AnalogInput1"].attrs[
        "long_name"
    ] = "External analog input 1 (fluorometer 1 in LISST-HAB & LISST-Black)"

    ds["LaserReferenceSensor"].attrs["units"] = "mW"
    ds["LaserReferenceSensor"].attrs["long_name"] = "Laser reference sensor"

    ds["Depth"].attrs["units"] = "m"
    ds["Depth"].attrs["standard_name"] = "depth"
    ds["Depth"].attrs["positive"] = "down"

    ds["Temperature"].attrs["units"] = "degree_C"
    ds["Temperature"].attrs["standard_name"] = "sea_water_temperature"

    ds["AnalogInput2"].attrs["units"] = "V"
    ds["AnalogInput2"].attrs[
        "long_name"
    ] = "External analog input 2 (fluorometer 2 in LISST-HAB & LISST-Black)"

    ds["MeanDiameter"].attrs["units"] = "micron"
    ds["MeanDiameter"].attrs[
        "long_name"
    ] = "Mean diameter (calculated from fully processed size distribution)"

    ds["TotalVolumeConcentration"].attrs["units"] = "ppm"
    ds["TotalVolumeConcentration"].attrs[
        "long_name"
    ] = "Total volume concentration (calculated from fully processed size distribution)"

    ds["RelativeHumidity"].attrs["units"] = "percent"
    ds["RelativeHumidity"].attrs["standard_name"] = "relative_humidity"
    ds["RelativeHumidity"].encoding["dtype"] = "i4"

    # ds["AccelerometerX"].attrs['units'] =
    ds["AccelerometerX"].attrs[
        "long_name"
    ] = "Accelerometer X [not presently calibrated or used]"
    ds["AccelerometerX"].encoding["dtype"] = "i4"

    # ds["AccelerometerY"].attrs['units'] =
    ds["AccelerometerY"].attrs[
        "long_name"
    ] = "Accelerometer Y [not presently calibrated or used]"
    ds["AccelerometerY"].encoding["dtype"] = "i4"

    # ds["AccelerometerZ"].attrs['units'] =
    ds["AccelerometerZ"].attrs[
        "long_name"
    ] = "Accelerometer Z [not presently calibrated or used]"
    ds["AccelerometerZ"].encoding["dtype"] = "i4"

    # ds["RawPressureMSB"].attrs['units'] =
    ds["RawPressureMSB"].attrs["long_name"] = "Raw pressure [most significant bit]"
    ds["RawPressureMSB"].encoding["dtype"] = "i4"

    # ds["RawPressureLSBs"].attrs['units'] =
    ds["RawPressureLSBs"].attrs[
        "long_name"
    ] = "Raw pressure [least significant 16 bits]"
    ds["RawPressureLSBs"].encoding["dtype"] = "i4"

    ds["AmbientLight"].attrs["units"] = "counts"
    ds["AmbientLight"].attrs["long_name"] = "Ambient light"
    ds["AmbientLight"].encoding["dtype"] = "i4"

    ds["AnalogInput3"].attrs["units"] = "V"
    ds["AnalogInput3"].attrs[
        "long_name"
    ] = "External analog input 3 (fluorometer 3 in LISST-HAB & LISST-Black)"

    ds["ComputedOpticalTransmissionOverPath"].attrs["units"] = "1"
    ds["ComputedOpticalTransmissionOverPath"].attrs[
        "long_name"
    ] = "Computed optical transmission over path"

    ds["BeamAttenuation"].attrs["units"] = "m-1"
    ds["BeamAttenuation"].attrs["long_name"] = "Beam attenuation (c)"

    ds["time"].attrs["standard_name"] = "time"

    ds["ring"].attrs["long_name"] = "Ring number"
    ds["ring"].encoding["dtype"] = "i4"

    ds["vc"].attrs["units"] = "uL L-1"
    ds["vc"].attrs["long_name"] = "Volume concentration"

    return ds


def get_ringsizes():
    """
    Return LISST ring sizes. Values obtained from:

    LISST-200X
    Particle Size Analyzer
    Including
    LISST-HAB & LISST-Black
    User’s Manual
    Version 2.3
    March, 2022

    Appendix B: Particle Size Bins (pp. 79)

    "There are 36 size ranges logarithmically placed from 1.00 – 500 microns in diameter.
    The upper size in each bin is approximately 1.18 times the lower, with the exception of bin 1.
    The table below shows the lower and upper limit of each size bin in microns, together with the median size (also in microns) for each size bin.
    The sizes are the same for both Spherical and Randomly Shaped inversions."
    """
    ringsizemedian = [
        1.21,
        1.60,
        1.89,
        2.23,
        2.63,
        3.11,
        3.67,
        4.33,
        5.11,
        6.03,
        7.11,
        8.39,
        9.90,
        11.7,
        13.8,
        16.3,
        19.2,
        22.7,
        26.7,
        31.6,
        37.2,
        43.9,
        51.9,
        61.2,
        72.2,
        85.2,
        101,
        119,
        140,
        165,
        195,
        230,
        273,
        324,
        386,
        459,
    ]
    ringsizeupper = [
        1.48,
        1.74,
        2.05,
        2.42,
        2.86,
        3.38,
        3.98,
        4.70,
        5.55,
        6.55,
        7.72,
        9.12,
        10.8,
        12.7,
        15.0,
        17.7,
        20.9,
        24.6,
        29.1,
        34.3,
        40.5,
        47.7,
        56.3,
        66.5,
        78.4,
        92.6,
        109,
        129,
        152,
        180,
        212,
        250,
        297,
        354,
        420,
        500,
    ]
    ringsizelower = [
        1.00,
        1.48,
        1.74,
        2.05,
        2.42,
        2.86,
        3.38,
        3.98,
        4.70,
        5.55,
        6.55,
        7.72,
        9.12,
        10.8,
        12.7,
        15.0,
        17.7,
        20.9,
        24.6,
        29.1,
        34.3,
        40.5,
        47.7,
        56.3,
        66.5,
        78.4,
        92.6,
        109,
        129,
        152,
        180,
        212,
        250,
        297,
        354,
        420,
    ]

    RSmedian = xr.DataArray(ringsizemedian, dims="ring")
    RSmedian.attrs["long_name"] = "particle size bin (median)"
    RSmedian.attrs["units"] = "micron"

    RSlower = xr.DataArray(ringsizelower, dims="ring")
    RSlower.attrs["long_name"] = "particle size bin (lower)"
    RSlower.attrs["units"] = "micron"

    RSupper = xr.DataArray(ringsizeupper, dims="ring")
    RSupper.attrs["long_name"] = "particle size bin (upper)"
    RSupper.attrs["units"] = "micron"

    return RSmedian, RSlower, RSupper


def check_and_reshape_burst(ds):
    if "operating_mode" in ds.attrs and ds.attrs["operating_mode"].lower() == "burst":
        dt0 = ds.time[1] - ds.time[0]
        spb = np.where(ds.time.diff(dim="time") != dt0)[0][0] + 1

        print(f"Inferred sample interval of {dt0.values} from timestamps")
        print(f"Inferred {spb} samples per burst from timestamps")

        r = np.shape(ds.time)[0]
        mod = r % spb

        if mod:
            print(
                "Number of rows is not a multiple of samples per burst; truncating to last full burst"
            )
            ds = ds.isel(time=range(r - mod))

        newtime = ds.time[0::spb]

        for v in ds.data_vars:
            if ds[v].dims == ("time",):
                attrsbak = ds[v].attrs
                ds[v] = xr.DataArray(
                    np.reshape(ds[v].values, (len(newtime), spb)),
                    dims=("newtime", "sample"),
                    coords=[
                        newtime,
                        range(spb),
                    ],
                )
                ds[v].attrs = attrsbak
            elif ds[v].dims == (
                "time",
                "ring",
            ):
                attrsbak = ds[v].attrs
                ds[v] = xr.DataArray(
                    np.reshape(ds[v].values, (len(newtime), spb, len(ds["ring"]))),
                    dims=("newtime", "sample", "ring"),
                    coords=[newtime, range(spb), ds["ring"]],
                )
                ds[v].attrs = attrsbak

        if "sample" in ds:
            ds["sample"].attrs["long_name"] = "Sample number"
            ds["sample"].attrs["units"] = "1"
            ds["sample"].encoding["dtype"] = "i4"

        ds = ds.drop("time")

        ds = ds.rename({"newtime": "time"})

    else:
        print("No operating_mode attr found; assuming fixed sample rate")

    return ds


def read_lop(ds):

    if "basefile_lop" not in ds.attrs:
        warnings.warn(
            "No .lop file specified; not including any metadata that would be in that file"
        )
        return ds

    with open(ds.attrs["basefile_lop"] + ".lop") as lop:
        for line in lop:
            ll = line.split(":", maxsplit=1)
            if len(ll) > 1:
                key = ll[0].replace(" ", "")
                value = ll[1].strip()
                ds.attrs[f"LISST{key}"] = value

    return ds
