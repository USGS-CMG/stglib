import numpy as np
import pandas as pd
import xarray as xr

from ..core import utils
from . import aqdutils


def hdr_to_cdf(metadata):
    """Load Aquadopp text files and output to netCDF format"""

    # TODO: clock drift code
    # TODO: logmeta code

    basefile = metadata["basefile"]

    if "prefix" in metadata:
        basefile = metadata["prefix"] + basefile

    utils.check_valid_globalatts_metadata(metadata)
    aqdutils.check_valid_config_metadata(metadata)

    # get instrument metadata from the HDR file
    instmeta = aqdutils.read_aqd_hdr(basefile)

    metadata["instmeta"] = instmeta

    print("Loading ASCII files")

    ds = xr.Dataset()

    # write out metadata first, then deal exclusively with xarray attrs
    ds = utils.write_metadata(ds, metadata)

    del metadata
    del instmeta

    # Load sensor data
    ds = load_sen(ds)

    ds = utils.ensure_cf(ds)

    # Deal with metadata peculiarities
    ds = aqdutils.check_attrs(ds, hr=True)

    ds = aqdutils.create_bindist(ds)

    # Load amplitude and velocity data
    ds = load_amp_vel_cor(ds, basefile)

    # Compute time stamps
    # ds = utils.shift_time(ds, ds.attrs["AQDAverageInterval"] / 2)

    # configure file
    if "prefix" in ds.attrs:
        cdf_filename = ds.attrs["prefix"] + ds.attrs["filename"] + "-raw.cdf"
    else:
        cdf_filename = ds.attrs["filename"] + "-raw.cdf"

    ds = aqdutils.update_attrs(ds, hr=True)

    # need to drop datetime
    # ds = ds.drop("datetime")

    ds.to_netcdf(cdf_filename, unlimited_dims=["time"])

    print("Finished writing data to %s" % cdf_filename)

    return ds


def load_sen(ds):
    """Load data from .sen file"""

    senfile = ds.attrs["basefile"] + ".sen"

    SEN = pd.read_csv(
        senfile,
        header=None,
        delim_whitespace=True,
        parse_dates={"datetime": [2, 0, 1, 3, 4, 5]},
        date_parser=aqdutils.date_parser,
        usecols=[0, 1, 2, 3, 4, 5, 6, 7, 10, 11, 12, 13, 14, 15, 16, 17, 18],
    )

    # rename columns from numeric to human-readable
    SEN.rename(
        columns={
            6: "Burst",
            7: "Ensemble",
            12: "Heading",
            13: "Pitch",
            14: "Roll",
            15: "Pressure",
            16: "Temperature",
            10: "Battery",
            11: "Soundspeed",
        },
        inplace=True,
    )

    SEN.rename(columns={17: "AnalogInput1"}, inplace=True)
    SEN["AnalogInput1"] = SEN["AnalogInput1"] * 5 / 65535
    SEN.rename(columns={18: "AnalogInput2"}, inplace=True)
    SEN["AnalogInput2"] = SEN["AnalogInput2"] * 5 / 65535

    SEN = SEN.rename(columns={"datetime": "time"}).set_index("time")

    RAW = SEN.to_xarray()

    r = np.shape(RAW.Heading)[0]
    mod = r % ds.attrs["AQDHRSamplesPerBurst"]
    if mod:
        print(
            "Number of rows is not a multiple of samples_per_burst; truncating to last full burst"
        )
        RAW = RAW.sel(time=RAW.time[0:-mod])

    ds["time"] = xr.DataArray(
        RAW.time[0 :: int(ds.attrs["AQDHRSamplesPerBurst"])].values, dims="time"
    )
    ds["time"] = pd.DatetimeIndex(ds["time"])

    ds["sample"] = xr.DataArray(range(ds.attrs["AQDHRSamplesPerBurst"]), dims="sample")

    for v in RAW.data_vars:
        ds[v] = xr.DataArray(
            np.reshape(RAW[v].values, (-1, int(ds.attrs["AQDHRSamplesPerBurst"]))),
            dims=["time", "sample"],
        )

    return ds


def load_amp_vel_cor(RAW, basefile):
    """Load amplitude, correlation, and velocity data from the .aN, .cN and .vN files"""

    for n in [1, 2, 3]:
        afile = basefile + ".a" + str(n)
        a = pd.read_csv(afile, header=None, delim_whitespace=True)

        spb = int(RAW.attrs["AQDHRSamplesPerBurst"])

        if "bindist" in RAW:
            coords = [RAW["time"], range(spb), RAW["bindist"]]
        else:
            coords = [RAW["time"], range(spb), RAW.attrs["AQDCCD"]]

        newshape = (len(coords[0]), len(coords[1]), len(coords[2]))

        RAW["AMP" + str(n)] = xr.DataArray(
            np.reshape(a.values[:, 2:], newshape),
            dims=("time", "sample", "bindist"),
            coords=coords,
        )

        cfile = basefile + ".c" + str(n)
        c = pd.read_csv(cfile, header=None, delim_whitespace=True)

        spb = int(RAW.attrs["AQDHRSamplesPerBurst"])

        if "bindist" in RAW:
            coords = [RAW["time"], range(spb), RAW["bindist"]]
        else:
            coords = [RAW["time"], range(spb), RAW.attrs["AQDCCD"]]

        newshape = (len(coords[0]), len(coords[1]), len(coords[2]))

        RAW["COR" + str(n)] = xr.DataArray(
            np.reshape(c.values[:, 2:], newshape),
            dims=("time", "sample", "bindist"),
            coords=coords,
        )

        vfile = basefile + ".v" + str(n)
        v = pd.read_csv(vfile, header=None, delim_whitespace=True)

        if RAW.attrs["AQDHRCoordinateSystem"] == "BEAM":
            thevars = {1: "VEL1", 2: "VEL2", 3: "VEL3"}
        elif RAW.attrs["AQDHRCoordinateSystem"] == "ENU":
            thevars = {1: "U", 2: "V", 3: "W"}
        elif RAW.attrs["AQDHRCoordinateSystem"] == "XYZ":
            thevars = {1: "X", 2: "Y", 3: "Z"}

        RAW[thevars[n]] = xr.DataArray(
            np.reshape(v.values[:, 2:], newshape),
            dims=("time", "sample", "bindist"),
            coords=coords,
        )

    return RAW
