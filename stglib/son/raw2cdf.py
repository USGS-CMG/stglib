import os
from math import nan

import numpy as np
import pandas as pd
import xarray as xr
from tqdm import tqdm

from ..core import utils
from . import sonutils


def read_81R(fname):
    """
    Parse an 81R file returned from an Imagenex Sonar 881A-GS
    """

    # Open and read the file in a binary format
    with open(fname, "rb") as fid:
        data = fid.read()

    # Universal lengths of the file header and the device list
    PINGHEADERBYTES = 1024
    DEVICELISTBYTES = 1024

    # Parse the start of the file
    time, header = sonutils.parse_pingHeader(data[:PINGHEADERBYTES], fname[-12:-8])

    # Gather the size of the switch commands and the return data header
    switchCommandBytes = header["SONSwitchCommandBytes"]
    returnHeaderBytes = header["SONReturnHeaderBytes"]

    # Drop unneeded variables
    header.pop("SONSwitchCommandBytes")
    header.pop("SONReturnHeaderBytes")
    header.pop("SONReturnDataHeaderType")

    # Extract information from the first set of switch commands
    offset = PINGHEADERBYTES + DEVICELISTBYTES
    SwitchSettings = sonutils.parse_switchCommand(
        data[offset : offset + switchCommandBytes]
    )

    # Add switch settings to header info
    header.update(SwitchSettings)

    # Calculate the number of pings using the file size and the size of each ping
    npings = int(os.path.getsize(fname) / header["SONTotalBytes"])

    # Convert the binary file into a 2-d numpy array
    imagedata = np.array(bytearray(data)).reshape(npings, header["SONTotalBytes"])

    # Initialize variables to make the loop more efficient
    variables = {}
    variables["ReturnDataHeaderType"] = [""] * npings
    # variables["HeadID"] = [""] * npings
    variables["HeadPosition"] = [nan] * npings
    variables["HeadAngle"] = [nan] * npings
    variables["StepDirection"] = [nan] * npings
    # variables["Range"] = [0] * npings
    variables["ProfileRange"] = [nan] * npings
    variables["NDataBytes"] = [nan] * npings
    variables["SonarPosition"] = [nan] * npings
    variables["SonarAngle"] = [nan] * npings
    variables["Pitch"] = [nan] * npings
    variables["Roll"] = [nan] * npings
    variables["Heading"] = [nan] * npings
    variables["NReturnBytes"] = [nan] * npings
    variables["GyroHeading"] = [nan] * npings

    # Extract data from each ping
    for i in range(npings):
        # Gather data from the universal header
        time, PingHeader = sonutils.parse_pingHeader(
            imagedata[i, :PINGHEADERBYTES].tobytes(), fname[-12:-8]
        )

        # Gather data from the switch commands
        offset = PINGHEADERBYTES + DEVICELISTBYTES
        SwitchCommand = sonutils.parse_switchCommand(
            imagedata[i, offset : offset + switchCommandBytes].tobytes()
        )

        # Gather data from the return data header
        offset += switchCommandBytes
        ReturnHeader = sonutils.parse_returnHeader(
            imagedata[i, offset : offset + returnHeaderBytes].tobytes(), SwitchCommand
        )

        # If the header isn't empty, transfer the gathered data into the output dictionary
        if ReturnHeader:
            variables["ReturnDataHeaderType"][i] = ReturnHeader["ReturnDataHeaderType"]
            # variables["HeadID"][i] = ReturnHeader["HeadID"]
            variables["HeadPosition"][i] = float(ReturnHeader["HeadPosition"])
            variables["HeadAngle"][i] = ReturnHeader["HeadAngle"]
            variables["StepDirection"][i] = float(ReturnHeader["StepDirection"])
            # variables["Range"][i] = ReturnHeader["Range"]
            variables["ProfileRange"][i] = float(ReturnHeader["ProfileRange"])
            variables["NDataBytes"][i] = ReturnHeader["NDataBytes"]
            variables["SonarPosition"][i] = ReturnHeader["SonarPosition"]
            variables["SonarAngle"][i] = ReturnHeader["SonarAngle"]
            variables["Pitch"][i] = ReturnHeader["Pitch"]
            variables["Roll"][i] = ReturnHeader["Roll"]
            variables["Heading"][i] = ReturnHeader["Heading"]
            variables["GyroHeading"][i] = ReturnHeader["GyroHeading"]
            # variables["time"] = time
            # variables['StartGain'][i] = SwitchCommand['StartGain']
            if "NReturnBytes" in ReturnHeader:
                variables["NReturnBytes"][i] = ReturnHeader["NReturnBytes"]
            else:
                variables["NReturnBytes"][i] = float("NaN")
                print(f"Problem at ping {i}")

    # Isolate the section of the numpy array that corresponds to the return data
    offset = PINGHEADERBYTES + DEVICELISTBYTES + switchCommandBytes + returnHeaderBytes
    image = imagedata[:, offset:-1]

    # Convert to xarray
    header.update({"SONReturnDataHeaderType": variables["ReturnDataHeaderType"][0]})

    variables.pop("ReturnDataHeaderType")
    variables.pop("NDataBytes")
    variables.pop("NReturnBytes")

    df = pd.DataFrame.from_dict(variables)
    df.index.names = ["scan"]
    # Reset index and add 1 to start scans from 1 instead of 0
    df.index = df.index.astype("int32") + 1
    ds = df.to_xarray()

    # Make xarray with time and echo data
    ds["points"] = np.array(range(1, image.shape[1] + 1), dtype="int32")
    echo_data = xr.DataArray(image, dims=["scan", "points"], name="sonar_image")
    ds = xr.merge([ds, echo_data])

    return ds, header, time


# Make raw CDF
def file81R_to_cdf(metadata):

    folder = metadata["basefile"]

    # List all files in the current directory
    files = os.listdir(folder)

    # Find unique number of sweeps
    sweepID = [file[6:8] for file in files]
    unique_sweep = list(set(sweepID))

    # Trim incomplete sweeps
    extra_files = len(files) % len(unique_sweep)
    if extra_files > 0:
        files = files[0:-extra_files]
        txt = "Trimmed incomplete sweeps"
        print(f"{txt}")

    # Find unique names for sets of sweeps
    names = [file[:-6] for file in files]
    unique_list = list(set(names))
    unique_list.sort()

    # For each sequence (5m or 20m)
    # Make progress bar with tqdm
    date = []
    sweep4_data = []
    for k in tqdm(
        range(0, len(unique_list)), bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt}"
    ):

        # Reset variables
        sweep_data = []
        ds_4sweeps = []

        # Find sets of sweeps
        sweep_list = [s for s in files if unique_list[k] in s]

        # Read header info from first sweep
        _, header, first_time = read_81R(folder + sweep_list[0])
        date.append(first_time)

        # Read each sweep in set
        for j in range(0, len(sweep_list)):
            sweep_data.append(read_81R(folder + sweep_list[j])[0])

        # Combine 4 sweeps
        ds_4sweeps = xr.concat(sweep_data, dim="sweep")

        # Add 4 sweeps to a list
        sweep4_data.append(ds_4sweeps)

    # Concatenate 4 sweep data
    ds = xr.concat(sweep4_data, dim="time")
    ds = ds.assign_coords(
        sweep=("sweep", np.array(range(1, len(unique_sweep) + 1), dtype="int32"))
    )
    ds = ds.assign_coords(time=("time", date))

    # Convert to int32
    for var in ds.data_vars:
        if ds[var].dtype == "int64" or ds[var].dtype == "uint8":
            ds[var] = ds[var].astype("int32")

    # Reorder dimensions
    ds = ds.transpose("time", "sweep", "scan", "points")

    ds = utils.write_metadata(ds, metadata)

    # Add trim sweep to history
    if extra_files > 0:
        ds = utils.insert_history(ds, txt)

    # Sort header alphabetically and add to global attributes
    header = sorted(header.items())
    ds.attrs.update(header)

    ds = utils.ensure_cf(ds)

    # Configure file
    cdf_filename = f"{ds.attrs['filename']}_{str(ds.attrs['SONRange'])}m-raw.cdf"

    ds.to_netcdf(cdf_filename, unlimited_dims=["time"], compute=False)

    print(f"Finished writing data to {cdf_filename}")

    return ds
