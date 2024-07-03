import warnings

import numpy as np
import pandas as pd
import xarray as xr

from ..aqd import aqdutils
from ..core import utils


def dat_to_cdf(metadata):
    basefile = metadata["basefile"]

    if "prefix" in metadata:
        basefile = metadata["prefix"] + basefile

    utils.check_valid_globalatts_metadata(metadata)
    aqdutils.check_valid_config_metadata(metadata, inst_type="VEC")

    # get instrument metadata from the HDR file
    instmeta = read_vec_hdr(basefile)

    metadata["instmeta"] = instmeta

    print("Loading ASCII files")

    dsvhd = load_vhd(basefile)

    ds = load_dat(basefile, metadata, dsvhd)

    ds = ds_checksum_check(ds)

    ds = utils.write_metadata(ds, metadata)

    ds = utils.ensure_cf(ds)

    ds = aqdutils.check_attrs(ds, inst_type="VEC")

    dssen = load_sen(basefile)

    # r = np.shape(dssen.Heading)[0]
    # senburstlen = int(ds.attrs["VECSamplesPerBurst"] / ds.attrs["VECSamplingRate"] + 1)
    ds.attrs["sample_interval"] = 1 / ds.attrs["VECSamplingRate"]

    # mod = r % senburstlen
    # if mod:
    # print(f"{ds.attrs['VECSamplesPerBurst']=}")
    # print(f"{ds.attrs['VECSamplingRate']=}")
    # histtext = "Number of rows in .sen file is not a multiple of ds.attrs['VECSamplesPerBurst']/ds.attrs['VECSamplingRate'] + 1; truncating to last full burst"
    # ds = utils.insert_history(ds, histtext)
    # dssen = dssen.sel(time=dssen.time[0:-mod])

    # dssen["timenew"] = xr.DataArray(dssen["time"].values[::senburstlen], dims="timenew")
    # dssen["sensample"] = xr.DataArray(range(senburstlen), dims="sensample")
    # for var in ["Heading", "Pitch", "Roll", "Battery", "Temperature", "Soundspeed", "ErrorCode", "StatusCode"]:
    #     dssen[var + "new"] = xr.DataArray(
    #         dssen[var].values.reshape((-1, senburstlen)),
    #         dims=["timenew", "sensample"],
    #     ).mean(dim="sensample")
    # for var in dssen.data_vars:
    #     if "new" not in var:
    #         dssen = dssen.drop(var)
    # for var in dssen.data_vars:
    #     dssen = dssen.rename({var: var.replace("new", "")})
    # dssen = dssen.drop(["time", "sensample"])
    # dssen = dssen.rename({"timenew": "time"})

    # Apply time from VHD file to DAT data
    # ds["time"] = dsvhd["time"]

    ds = ds.rename({"Ensemble": "sample"})

    # bring in sample volume distance from VHD, if burst mode
    if ds.attrs["VECBurstInterval"] != "CONTINUOUS":
        for var in [
            "DistProbeStartAvg",
            "DistSVolStartAvg",
            "DistProbeEndAvg",
            "DistSVolEndAvg",
        ]:
            ds[var] = dsvhd[var].reindex_like(ds, method="nearest")

    # ds = ds.swap_dims({"Burst": "time"})

    dssen = dssen.drop(
        ["Month", "Day", "Year", "Hour", "Minute", "Second", "AnalogInput"]
    )

    for var in dssen:
        ds[var] = dssen[var].reindex_like(ds, method="nearest")

    ds["TransMatrix"] = xr.DataArray(ds.attrs["VECTransMatrix"], dims=["inst", "beam"])
    ds["TransMatrix"].attrs["long_name"] = "Beam to XYZ (inst) Transformation Matrix"
    ds["TransMatrix"].attrs["note"] = "Provided by Nortek, based on transducer geometry"

    ds["inst"] = ["X", "Y", "Z"]
    ds["inst"].attrs["units"] = "1"
    ds["inst"].attrs["long_name"] = "Inst Reference Frame"

    ds["beam"] = [1, 2, 3]
    ds["beam"].attrs["units"] = "1"
    ds["beam"].attrs["long_name"] = "Beam Reference Frame"

    # Remove VECTransMatrix from attrs for netCDF compliance
    ds.attrs.pop("VECTransMatrix")

    for v in ds.data_vars:
        # need to do this or else a "coordinates" attribute with value of "Burst" hangs around
        ds[v].encoding["coordinates"] = None
    ds.encoding["coordinates"] = None

    # Compute time stamps. Only apply clock error or clock drift at this step since we still have time, sample dims
    ds = utils.shift_time(ds, 0)

    if "prefix" in ds.attrs:
        cdf_filename = ds.attrs["prefix"] + ds.attrs["filename"] + "-raw.cdf"
    else:
        cdf_filename = ds.attrs["filename"] + "-raw.cdf"

    ds.to_netcdf(cdf_filename, unlimited_dims=["time"])

    print(f"Finished writing data to {cdf_filename}")


def load_vhd(basefile):
    print(f"Loading {basefile}.vhd")
    names = [
        "Month",
        "Day",
        "Year",
        "Hour",
        "Minute",
        "Second",
        "Burst",
        "nsamp",
        "Noise1",
        "Noise2",
        "Noise3",
        "NoiseCor1",
        "NoiseCor2",
        "NoiseCor3",
        "DistProbeStart1",
        "DistProbeStart2",
        "DistProbeStart3",
        "DistProbeStartAvg",
        "DistSVolStartAvg",
        "DistProbeEnd1",
        "DistProbeEnd2",
        "DistProbeEnd3",
        "DistProbeEndAvg",
        "DistSVolEndAvg",
    ]
    vhd = pd.read_csv(f"{basefile}.vhd", sep=r"\s+", header=None, names=names)
    vhd["time"] = pd.to_datetime(
        vhd[["Year", "Month", "Day", "Hour", "Minute", "Second"]]
    )
    vhd = vhd.set_index("time")
    return vhd.to_xarray()


def load_sen(basefile):
    print(f"Loading {basefile}.sen")
    names = [
        "Month",
        "Day",
        "Year",
        "Hour",
        "Minute",
        "Second",
        "ErrorCode",
        "StatusCode",
        "Battery",
        "Soundspeed",
        "Heading",
        "Pitch",
        "Roll",
        "Temperature",
        "AnalogInput",
        "Checksum",
    ]
    sen = pd.read_csv(
        f"{basefile}.sen",
        sep=r"\s+",
        header=None,
        names=names,
        converters={"ErrorCode": str, "StatusCode": str},
    )
    sen["time"] = pd.to_datetime(
        sen[["Year", "Month", "Day", "Hour", "Minute", "Second"]]
    )

    x = sen["StatusCode"].values
    status = np.zeros(len(sen["time"]))

    for i in range(len(sen["StatusCode"])):
        status[i] = int(x[i], 2)

    sen["StatusCode"] = np.uint(status)

    sen["orientation"] = (sen["StatusCode"] & 1).astype(int)

    sen = sen.set_index("time")
    return sen.to_xarray()


def load_dat(basefile, metadata, dsvhd):
    print(f"Loading {basefile}.dat")
    names = [
        "Burst",
        "Ensemble",
        "VEL1",
        "VEL2",
        "VEL3",
        "AMP1",
        "AMP2",
        "AMP3",
        "SNR1",
        "SNR2",
        "SNR3",
        "COR1",
        "COR2",
        "COR3",
        "Pressure",
        "AnalogInput1",
        "AnalogInput2",
        "Checksum",
    ]

    dat = pd.read_csv(f"{basefile}.dat", header=None, sep=r"\s+", names=names)

    if metadata["instmeta"]["VECBurstInterval"] != "CONTINUOUS":
        array = np.zeros(np.shape(dat.VEL1))

        # Use burst number from vhd and dat files to match time
        for i in dsvhd["Burst"].values:
            index = np.where(dat.Burst == i)
            array[index] = dsvhd["time"][dsvhd["Burst"] == i]

        # Use ensemble number (sample number) and frequency so time increases correctly throughout the burst
        array = pd.to_datetime(array)
        frequency = 1000 / metadata["instmeta"]["VECSamplingRate"]
        delta = frequency * (dat.Ensemble.values - 1)
        dat["time"] = array + delta.astype("timedelta64[ms]")
        dat = dat.set_index(["time"])

    elif metadata["instmeta"]["VECBurstInterval"] == "CONTINUOUS":
        array = np.zeros(np.shape(dat.VEL1))
        array[:] = dsvhd["time"]
        # for continuous data, the ensemble values max out at 65535 and then repeat starting at 1
        # need to get this into continuous ensemble values to correctly create time
        # first get difference between ensemble numbers, in case there are any samples skipped
        sample = dat.Ensemble.diff()
        # since the first value will be NaN (can't take diff of the first value), need to make diff 1
        # since the ensemble maxes out at 65535, need to make step from 65535 back to 1 a diff value of 1
        index = np.where(sample == -65535)[0]
        index = np.append(index, [0])
        sample[index] = 1
        # next take cummlative sum so that ensemble/sample numbers are continous and will include possible sample jumps
        sample = (sample.cumsum()) - 1
        # Use ensemble number (sample number) and frequency so time increases correctly throughout the dataset
        array = pd.to_datetime(array)
        frequency = 1000 / metadata["instmeta"]["VECSamplingRate"]
        delta = frequency * (sample)
        dat["time"] = array + delta.astype("timedelta64[ms]")
        dat = dat.set_index(["time"])

    ds = dat.to_xarray()

    if np.any(ds["Checksum"].isnull()):
        warnings.warn(
            "Found NaN values in checksum; this indicates short bursts. Proceed with caution."
        )

    return ds


def read_vec_hdr(basefile):
    hdrFile = basefile + ".hdr"

    with open(hdrFile) as f:
        row = ""
        Instmeta = {}

        while "Hardware configuration" not in row:
            row = f.readline().rstrip()
            if "Sampling rate" in row:
                idx = row.find(" Hz")
                Instmeta["VECSamplingRate"] = int(row[38:idx])
            if "Nominal velocity range" in row:
                Instmeta["VECNominalVelocityRange"] = row[38:]
            if "Burst interval" in row:
                if "CONTINUOUS" not in row:
                    idx = row.find(" sec")
                    Instmeta["VECBurstInterval"] = int(row[38:idx])
                else:
                    Instmeta["VECBurstInterval"] = row[38:]
            if "Samples per burst" in row:
                if "N/A" not in row:
                    Instmeta["VECSamplesPerBurst"] = int(row[38:])
                else:
                    Instmeta["VECSamplesPerBurst"] = row[38:]
            if "Sampling volume" in row:
                Instmeta["VECSamplingVolume"] = row[38:]
            if "Measurement load" in row:
                idx = row.find(" %")
                Instmeta["VECMeasurementLoad"] = int(row[38:idx])
            if "Transmit length" in row:
                Instmeta["VECTransmitLength"] = row[38:]
            if "Receive length" in row:
                Instmeta["VECReceiveLength"] = row[38:]
            if "Output sync" in row:
                Instmeta["VECOutputSync"] = row[38:]
            if "Analog output" in row:
                Instmeta["VECAnalogOutput"] = row[38:]
            if "Analog input 1" in row:
                Instmeta["VECAnalogInput1"] = row[38:]
            if "Analog input 2" in row:
                Instmeta["VECAnalogInput2"] = row[38:]
            if "Power output" in row:
                Instmeta["VECAnalogInput2"] = row[38:]
            if "Output format" in row:
                Instmeta["VECOutputFormat"] = row[38:]
            if "Velocity scaling" in row:
                Instmeta["VECVelocityScaling"] = row[38:]
            if "IMU mode" in row:
                Instmeta["VECIMUMode"] = row[38:]
            if "IMU data type" in row:
                Instmeta["VECIMUDataType"] = row[38:]
            if "IMU Mag digital filter size" in row:
                Instmeta["VECIMUMagDigitalFilterSize"] = row[38:]
            if "IMU Gyro-Accel digital filter size" in row:
                Instmeta["VECIMUGyro/AccelDigitalFilterSize"] = row[38:]
            if "Powerlevel" in row:
                Instmeta["VECPowerLevel"] = row[38:]
            if "Coordinate system" in row:
                Instmeta["VECCoordinateSystem"] = row[38:]
            if "Sound speed" in row:
                Instmeta["VECSoundSpeed"] = row[38:]
            if "Salinity" in row:
                Instmeta["VECSalinity"] = row[38:]
            if "Distance between pings" in row:
                Instmeta["VECDistanceBetweenPings"] = row[38:]
            if "Number of beams" in row:
                Instmeta["VECDNumberOfBeams"] = int(row[38:])
            if "Software version" in row:
                Instmeta["VECSoftwareVersion"] = row[38:]
            if "Deployment name" in row:
                Instmeta["VECDeploymentName"] = row[38:]
            if "Wrap mode" in row:
                Instmeta["VECWrapMode"] = row[38:]
            if "Deployment time" in row:
                Instmeta["VECDeploymentTime"] = row[38:]
            if "Comments" in row:
                Instmeta["VECComments"] = row[38:]
                # There may be up to three lines of comments, but only if they were added during deployment.
                # These extra lines will be preceded by blanks instead of a field name.
                # After the comments lines are the System lines, which we currenty don't handle,
                # so we can read the next lines and add to Comments if present.
                # Example showing a two-line comment followed by System1 field:
                # Comments                              SP 15916 30 cmab
                #                                       WTS21
                # System1                               19
                for n in range(2):
                    row = f.readline().rstrip()
                    if len(row) and row[0] == " ":
                        Instmeta["VECComments"] += "\n"
                        Instmeta["VECComments"] += row[38:]

        while "Head configuration" not in row:
            row = f.readline().rstrip()
            if "Serial number" in row:
                Instmeta["VECSerialNumber"] = row[38:]
            elif "Hardware revision" in row:
                Instmeta["VECHardwareRevision"] = row[38:]
            elif "Recorder size" in row:
                Instmeta["VECRecorderSize"] = row[38:]
            elif "Firmware version" in row:
                Instmeta["VECFirmwareVersion"] = row[38:]
            elif "Power output" in row:
                Instmeta["VECAnalogPowerOutput"] = row[38:]

        while "Data file format" not in row:
            row = f.readline().rstrip()
            if "Pressure sensor" in row:
                Instmeta["VECPressureSensor"] = row[38:]
            elif "Compass" in row:
                Instmeta["VECCompass"] = row[38:]
            elif "Tilt sensor" in row:
                Instmeta["VECTiltSensor"] = row[38:]
            elif "Head frequency" in row:
                idx = row.find(" kHz")
                Instmeta["VECFrequency"] = int(row[38:idx])
            elif "Serial number" in row:
                Instmeta["VECHeadSerialNumber"] = row[38:]
            elif "IMU Sensor" in row:
                Instmeta["VECIMUSensor"] = row[38:]
            elif "Transformation matrix" in row:
                Instmeta["VECTransMatrix"] = np.zeros((3, 3))
                Instmeta["VECTransMatrix"][0, :] = [float(x) for x in row[38:].split()]
                row = f.readline().rstrip()
                Instmeta["VECTransMatrix"][1, :] = [float(x) for x in row[38:].split()]
                row = f.readline().rstrip()
                Instmeta["VECTransMatrix"][2, :] = [float(x) for x in row[38:].split()]
            elif "Pressure sensor calibration" in row:
                Instmeta["VECPressureCal"] = row[38:]

    return Instmeta


def ds_checksum_check(ds):
    if np.any(ds["Checksum"] == 1):
        warnings.warn(
            "Non-zero checksum values found in data. This indicates a failed checksum and potentially bad data. Proceed with caution."
        )

    return ds
