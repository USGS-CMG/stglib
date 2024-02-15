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
    aqdutils.check_valid_config_metadata(metadata)

    # get instrument metadata from the HDR file
    instmeta = read_vec_hdr(basefile)

    metadata["instmeta"] = instmeta

    print("Loading ASCII files")

    ds = load_dat(basefile)

    ds = ds_checksum_check(ds)

    ds = utils.write_metadata(ds, metadata)

    ds = utils.ensure_cf(ds)

    ds = aqdutils.check_attrs(ds, inst_type="VEC")

    dsvhd = load_vhd(basefile)

    dssen = load_sen(basefile)

    r = np.shape(dssen.Heading)[0]
    senburstlen = int(ds.attrs["VECSamplesPerBurst"] / ds.attrs["VECSamplingRate"] + 1)
    ds.attrs["sample_interval"] = 1 / ds.attrs["VECSamplingRate"]

    mod = r % senburstlen
    if mod:
        print(f"{ds.attrs['VECSamplesPerBurst']=}")
        print(f"{ds.attrs['VECSamplingRate']=}")
        print(
            "Number of rows in .sen file is not a multiple of ds.attrs['VECSamplesPerBurst']/ds.attrs['VECSamplingRate'] + 1; truncating to last full burst"
        )
        dssen = dssen.sel(time=dssen.time[0:-mod])

    dssen["timenew"] = xr.DataArray(dssen["time"].values[::senburstlen], dims="timenew")
    dssen["sensample"] = xr.DataArray(range(senburstlen), dims="sensample")
    for var in ["Heading", "Pitch", "Roll", "Battery", "Temperature", "Soundspeed"]:
        dssen[var + "new"] = xr.DataArray(
            dssen[var].values.reshape((-1, senburstlen)),
            dims=["timenew", "sensample"],
        ).mean(dim="sensample")
    for var in dssen.data_vars:
        if "new" not in var:
            dssen = dssen.drop(var)
    for var in dssen.data_vars:
        dssen = dssen.rename({var: var.replace("new", "")})
    dssen = dssen.drop(["time", "sensample"])
    dssen = dssen.rename({"timenew": "time"})

    # Apply time from VHD file to DAT data
    ds["time"] = dsvhd["time"]

    ds = ds.swap_dims({"Burst": "time"})
    ds = ds.rename({"Ensemble": "sample"})

    for var in ["Heading", "Pitch", "Roll", "Battery", "Temperature", "Soundspeed"]:
        ds[var] = dssen[var].reindex_like(ds, method="nearest")

    ds["TransMatrix"] = xr.DataArray(ds.attrs["VECTransMatrix"])
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
    vhd = pd.read_csv(
        f"{basefile}.vhd", delim_whitespace=True, header=None, names=names
    )
    vhd["time"] = pd.to_datetime(
        vhd[["Year", "Month", "Day", "Hour", "Minute", "Second"]]
    )
    vhd = vhd.set_index("Burst")
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
        f"{basefile}.sen", delim_whitespace=True, header=None, names=names
    )
    sen["time"] = pd.to_datetime(
        sen[["Year", "Month", "Day", "Hour", "Minute", "Second"]]
    )
    sen = sen.set_index("time")
    return sen.to_xarray()


def load_dat(basefile):
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
    dat = pd.read_csv(
        f"{basefile}.dat", header=None, delim_whitespace=True, names=names
    )
    dat = dat.set_index(["Burst", "Ensemble"])
    return dat.to_xarray()


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
                idx = row.find(" sec")
                Instmeta["VECBurstInterval"] = int(row[38:idx])
            if "Samples per burst" in row:
                Instmeta["VECSamplesPerBurst"] = int(row[38:])
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
