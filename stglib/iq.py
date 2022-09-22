import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from . import core
from .core import utils


def mat_to_cdf(metadata):
    """
    Process SonTek IQ .mat data to raw .cdf file
    """

    basefile = metadata["basefile"]

    ds = read_iq(basefile + ".mat")

    # write out metadata first, then deal exclusively with xarray attrs
    ds = utils.write_metadata(ds, metadata)

    del metadata

    # configure file
    cdf_filename = ds.attrs["filename"] + "-raw.cdf"

    ds.to_netcdf(cdf_filename, unlimited_dims=["time"])

    print("Finished writing data to %s" % cdf_filename)

    return ds


def read_iq(filnam):
    """Read SonTek IQ data which has been exported as a Matlab .mat file from IQ
    software into an xarray Dataset

    Parameters
    ----------
    filnam : string
        The SonTek .mat filename

    Returns
    -------
    xarray.Dataset
        An xarray Dataset of the IQ data
    """

    iqmat = core.utils.loadmat(filnam)
    # offset = iqmat['FlowSubData_PrfHeader_0_BlankingDistance']
    # beamdist_0 = np.linspace(offset, offset + \
    # 100*iqmat['FlowSubData_PrfHeader_0_CellSize'], 100)
    ds = {}

    ds["time"] = xr.DataArray(
        iqmat["FlowData_SampleTime"],
        attrs={
            "standard_name": "time",
            "axis": "T",
            # per email from SonTek
            "units": "microseconds since 2000-01-01 00:00:00",
            "calendar": "proleptic_gregorian",
        },
        dims="time",
    )

    ds["velbeam"] = xr.DataArray(
        [1, 2, 3, 4],
        dims="velbeam",
        attrs={"long_name": "velocity beam number", "units": "1"},
    )
    ds["beam"] = xr.DataArray(
        [1, 2, 3, 4, 5],
        dims="beam",
        attrs={"long_name": "beam number", "units": "1"},
    )
    # ds['beamdist_0'] = xr.DataArray(beamdist_0, dims='beamdist_0')
    # attrs = {}

    # need to do this because sometimes the flowsubdata and profile data is
    # one burst longer
    timelen = len(ds["time"])

    for k in iqmat:
        if "__" not in k and "FlowSubData" not in k:
            # print(k, np.shape(iqmat[k]))
            if len(np.ravel(iqmat[k])) == len(ds["time"]):
                ds[k] = xr.DataArray(np.ravel(iqmat[k]), dims="time")
                if k in iqmat["Data_Units"]:
                    ds[k].attrs["units"] = iqmat["Data_Units"][k].replace("/s", " s-1")
            elif "_2_" in k or "_3_" in k:
                ds[k] = xr.DataArray(
                    iqmat[k][0:timelen, :], dims=("time", "cell_across")
                )
                if k in iqmat["Data_Units"]:
                    ds[k].attrs["units"] = iqmat["Data_Units"][k].replace("/s", " s-1")
            elif "_0_" in k or "_1_" in k:
                ds[k] = xr.DataArray(
                    iqmat[k][0:timelen, :], dims=("time", "cell_along")
                )
                if k in iqmat["Data_Units"]:
                    ds[k].attrs["units"] = iqmat["Data_Units"][k].replace("/s", " s-1")
            elif "FlowData_Vel" in k or "FlowData_SNR" in k:
                ds[k] = xr.DataArray(iqmat[k][0:timelen, :], dims=("time", "velbeam"))
                if k in iqmat["Data_Units"]:
                    ds[k].attrs["units"] = iqmat["Data_Units"][k].replace("/s", " s-1")
            elif "FlowData_NoiseLevel" in k:
                ds[k] = xr.DataArray(iqmat[k][0:timelen, :], dims=("time", "beam"))
                if k in iqmat["Data_Units"]:
                    ds[k].attrs["units"] = iqmat["Data_Units"][k].replace("/s", " s-1")

    ds["cell_along"] = np.arange(ds["Profile_0_Vel"].shape[1])
    ds["cell_across"] = np.arange(ds["Profile_2_Vel"].shape[1])

    ds = xr.Dataset(ds)
    for k in iqmat["System_IqSetup"]["basicSetup"]:
        if "spare" not in k:
            ds.attrs[k] = iqmat["System_IqSetup"]["basicSetup"][k]
    for k in iqmat["System_Id"]:
        ds.attrs[k] = iqmat["System_Id"][k]
    for k in iqmat["System_IqState"]:
        if "spare" not in k:
            ds.attrs[k.replace("[", "_").replace("]", "_")] = iqmat["System_IqState"][k]

    return xr.decode_cf(ds)


def rename_vars(ds):

    # set up dict of instrument -> EPIC variable names
    varnames = {
        "Batt": "Bat_106",
        "Temp": "T_28",
        "Pitch": "Ptch_1216",
        "Roll": "Roll_1217",
        "Depth": "D_3",
        "Pressure": "P_1",
        "AdjustedPressure": "P_1ac",
    }

    # check to make sure they exist before trying to rename
    newvars = {}
    for k in varnames:
        if k in ds:
            newvars[k] = varnames[k]

    return ds.rename(newvars)


def remove_FlowData(ds):

    newvars = {}
    for k in ds:
        newvars[k] = k.replace("FlowData_", "")

    return ds.rename(newvars)


def clean_iq(iq):
    """
    Preliminary data cleaning when SNR < 0
    """

    iq["Vel_Mean"].values[iq["Vel_Mean"] < -214748] = np.nan
    iq["Vel"].values[iq["Vel"] == -214748368] = np.nan
    for bm in range(4):
        pr = "Profile_" + str(bm) + "_Vel"
        iq[pr].values[iq[pr] == -214748368] = np.nan
        am = "Profile_" + str(bm) + "_Amp"
        iq[am].values[iq[am] == 65535] = np.nan
        st = "Profile_" + str(bm) + "_VelStd"
        iq[st].values[iq[st] < 0] = np.nan

    return iq


def vel_to_ms(iq):
    """
    Convert velocity data from mm/s to m/s
    """

    for var in ["FlowData_Vel_Mean", "FlowData_Vel"]:
        iq[var] = iq[var] / 1000

    return iq


def make_beamdist(iq):
    """
    Generate physical coordinates to pair with the logical beamdist coordinates
    """
    for bm in range(4):
        if bm < 2:
            bdname = "beamdist_0_1"
        else:
            bdname = "beamdist_2_3"

        r = range(len(iq[bdname]))

        time = np.tile(iq["time"], (len(iq[bdname]), 1)).transpose()

        cells = np.zeros(np.shape(iq["Profile_" + str(bm) + "_Vel"]))
        fsdbd = "FlowSubData_PrfHeader_" + str(bm) + "_BlankingDistance"
        fsdcs = "FlowSubData_PrfHeader_" + str(bm) + "_CellSize"
        for n in range(len(iq["time"])):
            cells[n, :] = iq[fsdbd][n].values + r * iq[fsdcs][n].values
        iq["cells_" + str(bm)] = xr.DataArray(cells, dims=("time", bdname))
        iq["time_" + str(bm)] = xr.DataArray(time, dims=("time", bdname))
        iq = iq.set_coords(["cells_" + str(bm), "time_" + str(bm)])

    return iq


def make_iq_plots(iq, directory="", savefig=False):
    """
    Make IQ turnaround plots
    """

    plt.figure(figsize=(11, 8.5))

    for n, var in enumerate(
        ["FlowData_Depth", "FlowData_Vel_Mean", "FlowData_Flow"], start=1
    ):
        plt.subplot(3, 1, n)
        plt.plot(iq["time"], iq[var])
        plt.ylabel(var + " [" + iq[var].attrs["units"] + "]")

    if savefig:
        plt.savefig(directory + "/iq_stage_vel_flow.pdf")
    plt.show()


def cdf_to_nc(cdf_filename):
    """
    Load a "raw" .cdf file and generate a processed .nc file
    """

    # Load raw .cdf data
    ds = xr.open_dataset(cdf_filename)

    ds = remove_FlowData(ds)

    # Clip data to in/out water times or via good_ens
    ds = utils.clip_ds(ds)

    ds = clean_iq(ds)

    # assign min/max:
    ds = utils.add_min_max(ds)

    ds = utils.add_start_stop_time(ds)

    ds = utils.add_delta_t(ds)

    # add lat/lon coordinates
    ds = utils.ds_add_lat_lon(ds)

    ds = rename_vars(ds)

    ds = ds_add_attrs(ds)

    ds = utils.no_p_create_depth(ds)

    ds = ds.drop(["SampleNumber", "SampleTime"])

    # add lat/lon coordinates to each variable
    for var in ds.variables:
        if (var not in ds.coords) and ("time" not in var):
            # ds = utils.add_lat_lon(ds, var)
            # cast as float32
            ds = utils.set_var_dtype(ds, var)

    dsflow = ds.copy()
    dsprof = ds.copy()

    dsflow = dsflow.drop([k for k in dsflow if "Profile_" in k])
    dsprof = dsprof.drop([k for k in dsprof if "Profile_" not in k])

    # Write to .nc file
    print("Writing cleaned/trimmed data to .nc file")

    nc_filename = dsflow.attrs["filename"] + "flow-a.nc"
    dsflow.to_netcdf(nc_filename, format=format, unlimited_dims=["time"])
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])
    print("Done writing netCDF file", nc_filename)

    nc_filename = dsprof.attrs["filename"] + "prof-a.nc"
    dsprof.to_netcdf(nc_filename, format=format)
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])
    print("Done writing netCDF file", nc_filename)


def ds_add_attrs(ds):

    ds.attrs["serial_number"] = ds.attrs["SerialNumber"]
    ds.attrs["INST_TYPE"] = "SonTek-IQ Plus"

    # Update attributes for EPIC and STG compliance
    ds = utils.ds_coord_no_fillvalue(ds)

    ds["time"].attrs.update(
        {"standard_name": "time", "axis": "T", "long_name": "time (UTC)"}
    )

    if "epic_time" in ds:
        ds["epic_time"].attrs.update(
            {"units": "True Julian Day", "type": "EVEN", "epic_code": 624}
        )

    if "epic_time2" in ds:
        ds["epic_time2"].attrs.update(
            {"units": "msec since 0:00 GMT", "type": "EVEN", "epic_code": 624}
        )

    ds["D_3"].attrs.update(
        {"long_name": "Depth (relative to the top of the instrument)", "epic_code": 3}
    )

    # descriptions from Sontek-IQ Series User's Manual available at
    # http://info.xylem.com/sontek-iq-manual.html
    ds["Stage"].attrs["long_name"] = "Stage (water depth of the user-defined channel)"
    ds["Area"].attrs["long_name"] = "Cross-sectional area of user-defined channel"
    ds["Flow"].attrs["long_name"] = "Flow rate (using defined channel geometry)"

    ds["Vel_Mean"].attrs["long_name"] = "Mean velocity"
    ds["Volume_Total"].attrs[
        "long_name"
    ] = "Total water volume (based on all measured flow)"
    ds["Volume_Positive"].attrs[
        "long_name"
    ] = "Total volume of water in the positive downstream direction"
    ds["Volume_Negative"].attrs[
        "long_name"
    ] = "Total volume of water in the negative upstream direction"
    ds["Vel"].attrs["long_name"] = "Velocity"
    ds["VelXYZ"].attrs["long_name"] = ""
    ds["VelStd"].attrs["long_name"] = "Velocity standard deviation"
    ds["SNR"].attrs["long_name"] = "Signal-to-noise ratio"
    ds["NoiseLevel"].attrs["long_name"] = "Acoustic noise level"
    ds["Range"].attrs["long_name"] = "Acoustically measured distance to water surface"
    ds["T_28"].attrs["long_name"] = "Water temperature"
    ds["P_1"].attrs["long_name"] = "Pressure"
    ds["PressOffsetAdjust"].attrs[
        "long_name"
    ] = "Atmospheric pressure adjustment (see SonTek-IQ User's Manual for details)"
    ds["P_1ac"].attrs[
        "long_name"
    ] = "Measurement with atmospheric pressure removed (see SonTek-IQ User's Manual for details)"
    ds["Bat_106"].attrs["long_name"] = "Battery voltage"
    ds["Ptch_1216"].attrs["long_name"] = "Pitch angle in degrees"
    # to be UDUNITS compatible
    if ds["Ptch_1216"].attrs["units"] == "deg":
        ds["Ptch_1216"].attrs["units"] = "degree"
    if ds["Roll_1217"].attrs["units"] == "deg":
        ds["Roll_1217"].attrs["units"] = "degree"
    ds["Roll_1217"].attrs["long_name"] = "Roll angle in degrees"
    ds["VbPercentGood"].attrs["long_name"] = "Vertical beam percent good"
    ds["HorizontalSkew"].attrs["long_name"] = "Horizontal skew"
    ds["SystemInWater"].attrs[
        "long_name"
    ] = "Percentage of sample during which instrument was submerged (100% means it was submerged for entire sample)"

    # Profile Variables
    for n in range(4):
        ds["Profile_%d_Amp" % n].attrs["long_name"] = "Beam %d amplitude" % n
        ds["Profile_%d_VelStd" % n].attrs["long_name"] = (
            "Beam %d velocity profile standard deviation" % n
        )
        ds["Profile_%d_Vel" % n].attrs["long_name"] = "Beam %d velocity profile" % n

    def add_attributes(var, dsattrs):
        var.attrs.update(
            {
                "initial_instrument_height": dsattrs["initial_instrument_height"],
                # 'nominal_instrument_depth': dsattrs['nominal_instrument_depth'],
                "height_depth_units": "m",
            }
        )
        # var.encoding["_FillValue"] = 1e35

    for var in ds.variables:
        if (var not in ds.coords) and ("time" not in var):
            add_attributes(ds[var], ds.attrs)

    ds.attrs["COMPOSITE"] = np.int32(0)

    return ds
