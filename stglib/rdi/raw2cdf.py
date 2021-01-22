import pandas as pd
import xarray as xr
import matplotlib.dates

from ..core import utils
from . import rdradcp, rdiadcpy


def raw_to_cdf(metadata):
    """Load a Aquadopp text files and output to netCDF format"""

    # TODO: clock drift code
    # TODO: logmeta code

    basefile = metadata["basefile"]

    if "prefix" in metadata:
        basefile = metadata["prefix"] + basefile

    ds = xr.Dataset()

    # write out metadata first, then deal exclusively with xarray attrs
    ds = utils.write_metadata(ds, metadata)

    ensemble_count, netcdf_index, ens_error, alldata = rdiadcpy.convert_pd0_to_netcdf(
        basefile, [0, -1], "0", "CF", None
    )

    ds["time"] = xr.DataArray(pd.DatetimeIndex(alldata["time"]), dims="time")
    ds["bindist"] = xr.DataArray(alldata["bindist"], dims="bindist")
    for k in alldata:
        if k is "time" or k is "bindist" or k is "FLeader":
            continue
        if alldata[k].ndim == 1 and alldata[k] != []:
            ds[k] = xr.DataArray(alldata[k], dims="time")
        if alldata[k].ndim == 2:
            ds[k] = xr.DataArray(alldata[k], dims=["time", "bindist"])

    # varobj = cdf.createVariable('bindist', 'f4', ('depth',), fill_value=floatfill)
    # note name is one of the netcdf4 reserved attributes, use setncattr
    # varobj.setncattr('name', "bindist")
    ds["bindist"].attrs["units"] = "m"
    ds["bindist"].attrs["long_name"] = "bin distance from instrument for slant beams"
    ds["bindist"].attrs["epic_code"] = 0
    ds["bindist"].attrs[
        "NOTE"
    ] = "distance is calculated from center of bin 1 and bin size"

    # varobj = cdf.createVariable('Rec', 'u4', 'time', fill_value=intfill)
    ds["Rec"].attrs["units"] = "count"
    ds["Rec"].attrs["long_name"] = "Ensemble Number"
    # the ensemble number is a two byte LSB and a one byte MSB (for the rollover)
    # varobj.valid_range = [0, 2**23]

    for i in range(4):
        varname = f"vel{i+1}"
        ds[varname].attrs["units"] = "mm s-1"
        ds[varname].attrs["long_name"] = f"Beam {i+1} velocity"
        ds[varname].attrs["epic_code"] = 1277 + i
        ds[varname].encoding["dtype"] = 'i2'

        varname = "cor%d" % (i + 1)
        ds[varname].attrs["units"] = "counts"
        ds[varname].attrs["long_name"] = "Beam %d correlation" % (i + 1)
        ds[varname].attrs["epic_code"] = 1285 + i
        ds[varname].encoding['dtype'] = 'u2'

        varname = "att%d" % (i + 1)
        ds[varname].attrs["units"] = "counts"
        ds[varname].attrs["epic_code"] = 1281 + i
        ds[varname].attrs["long_name"] = "ADCP attenuation of beam %d" % (i + 1)
        ds[varname].encoding['dtype'] = 'u2'

        varname = "PGd%d" % (i + 1)
        if varname in ds:
            ds[varname].attrs["units"] = "counts"
            ds[varname].attrs["long_name"] = "Percent Good Beam %d" % (i + 1)
            ds[varname].attrs["epic_code"] = 1241 + i

        varname = "EWD%d" % (i + 1)
        ds[varname].attrs["units"] = "binary flag"
        ds[varname].attrs["long_name"] = "Error Status Word %d" % (i + 1)
        ds[varname].encoding["dtype"] = 'u2'

    # varobj = cdf.createVariable('Hdg', 'f4', ('time',), fill_value=floatfill)
    ds["Hdg"].attrs["units"] = "hundredths of degrees"
    ds["Hdg"].attrs["standard_name"] = "platform_orientation"
    ds["Hdg"].attrs["epic_code"] = 1215
    ds["Hdg"].attrs["heading_alignment"] = alldata["FLeader"][
        "Heading_Alignment_Hundredths_of_Deg"
    ]
    ds["Hdg"].attrs["heading_bias"] = alldata["FLeader"][
        "Heading_Bias_Hundredths_of_Deg"
    ]
    if alldata["FLeader"]["Heading_Bias_Hundredths_of_Deg"] == 0:
        ds["Hdg"].attrs[
            "NOTE_9"
        ] = "no heading bias was applied by EB during deployment or by wavesmon"
    else:
        ds["Hdg"].attrs[
            "NOTE_9"
        ] = "a heading bias was applied by EB during deployment or by wavesmon"

    # varobj = cdf.createVariable('Ptch', 'f4', ('time',), fill_value=floatfill)
    ds["Ptch"].attrs["units"] = "hundredths of degrees"
    ds["Ptch"].attrs["standard_name"] = "platform_pitch"
    ds["Ptch"].attrs["epic_code"] = 1216

    # varobj = cdf.createVariable('Roll', 'f4', ('time',), fill_value=floatfill)
    ds["Roll"].attrs["units"] = "hundredths of degrees"
    ds["Roll"].attrs["standard_name"] = "platform_roll"
    ds["Roll"].attrs["epic_code"] = 1217

    # varobj = cdf.createVariable('sv', 'f4', ('time',), fill_value=floatfill)
    ds["sv"].attrs["units"] = "m s-1"
    ds["sv"].attrs["standard_name"] = "speed_of_sound_in_sea_water"
    ds["sv"].encoding["dtype"] = 'u2'

    # varobj = cdf.createVariable('HdgSTD', 'f4', ('time',), fill_value=floatfill)
    ds["HdgSTD"].attrs["units"] = "degrees"
    ds["HdgSTD"].attrs["long_name"] = "Heading Standard Deviation"

    # varobj = cdf.createVariable('PtchSTD', 'f4', ('time',), fill_value=floatfill)
    ds["PtchSTD"].attrs["units"] = "tenths of degrees"
    ds["PtchSTD"].attrs["long_name"] = "Pitch Standard Deviation"

    # varobj = cdf.createVariable('RollSTD', 'f4', ('time',), fill_value=floatfill)
    ds["RollSTD"].attrs["units"] = "tenths of degrees"
    ds["RollSTD"].attrs["long_name"] = "Roll Standard Deviation"

    # varobj = cdf.createVariable('Tx', 'f4', ('time',), fill_value=floatfill)
    ds["Tx"].attrs["units"] = "hundredths of degrees"
    ds["Tx"].attrs["long_name"] = "ADCP Transducer Temperature"
    ds["Tx"].attrs["epic_code"] = 3017

    # varobj = cdf.createVariable('S', 'f4', ('time',), fill_value=floatfill)
    ds["S"].attrs["units"] = "PPT"
    ds["S"].attrs["standard_name"] = "sea_water_salinity"
    ds["S"].attrs["epic_code"] = 40

    # varobj'] = cdf.createVariable('xmitc', 'f4', ('time',), fill_value=floatfill)
    ds["xmitc"].attrs["units"] = "amps"
    ds["xmitc"].attrs["long_name"] = "transmit current"

    # varobj'] = cdf.createVariable('xmitv', 'f4', ('time',), fill_value=floatfill)
    ds["xmitv"].attrs["units"] = "volts"
    ds["xmitv"].attrs["long_name"] = "transmit voltage"

    # varobj'] = cdf.createVariable('Ambient_Temp', 'i2', ('time',), fill_value=intfill)
    ds["Ambient_Temp"].attrs["units"] = "C"
    ds["Ambient_Temp"].attrs["long_name"] = "Ambient_Temp"

    # varobj'] = cdf.createVariable('Pressure+', 'i2', ('time',), fill_value=intfill)
    ds["Pressure+"].attrs["units"] = "unknown"
    ds["Pressure+"].attrs["long_name"] = "Pressure+"

    # varobj'] = cdf.createVariable('Pressure-', 'i2', ('time',), fill_value=intfill)
    ds["Pressure-"].attrs["units"] = "unknown"
    ds["Pressure-"].attrs["long_name"] = "Pressure-"

    # varobj'] = cdf.createVariable('Attitude_Temp', 'i2', ('time',), fill_value=intfill)
    ds["Attitude_Temp"].attrs["units"] = "C"
    ds["Attitude_Temp"].attrs["long_name"] = "Attitude_Temp"

    if "Pressure" in ds:
        # varobj = cdf.createVariable('Pressure', 'f4', ('time',), fill_value=floatfill)
        ds["Pressure"].attrs["units"] = "deca-pascals"
        ds["Pressure"].attrs["long_name"] = "ADCP Transducer Pressure"
        ds["Pressure"].attrs["epic_code"] = 4

    if "PressVar" in ds:
        # varobj = cdf.createVariable('PressVar', 'f4', ('time',), fill_value=floatfill)
        ds["PressVar"].attrs["units"] = "deca-pascals"
        ds["PressVar"].attrs["long_name"] = "ADCP Transducer Pressure Variance"

    if "vel5" in ds:  # 5-beam instrument
        # varobj = cdf.createVariable("vel5", 'f4', ('time', 'depth'), fill_value=floatfill)
        ds["vel5"].attrs["units"] = "mm s-1"
        ds["vel5"].attrs["long_name"] = "Beam 5 velocity (mm s-1)"
        ds["vel5"].encoding["dtype"] = 'i2'

        ds["cor5"].attrs["units"] = "counts"
        ds["cor5"].attrs["long_name"] = "Beam 5 correlation"
        ds["cor5"].encoding["dtype"] = 'u2'

        ds["att5"].attrs["units"] = "counts"
        ds["att5"].attrs["long_name"] = "ADCP attenuation of beam 5"
        ds["att5"].encoding["dtype"] = 'u2'

    # add FLeader attrs
    for k in alldata["FLeader"]:
        if k not in ds.attrs:
            # print(f'adding {k}')
            ds.attrs[f"RDI_{k}"] = alldata["FLeader"][k]
        else:
            print(f'*** warning {k} already in global attributes ***')

    # configure file
    if "prefix" in ds.attrs:
        cdf_filename = ds.attrs["prefix"] + ds.attrs["filename"] + "-raw.cdf"
    else:
        cdf_filename = ds.attrs["filename"] + "-raw.cdf"

    ds.to_netcdf(cdf_filename, unlimited_dims=["time"])

    print("Finished writing data to %s" % cdf_filename)

    return ds
