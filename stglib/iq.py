import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from . import core
from .core import utils
from .core import qaqc


def mat_to_cdf(metadata):
    """
    Process SonTek IQ .mat data to raw .cdf file
    """

    basefile = metadata["basefile"]

    ds = read_iq(basefile + ".mat")
    
    ds = create_iqbindist(ds)

    # write out metadata first, then deal exclusively with xarray attrs
    ds = utils.write_metadata(ds, metadata)

    del metadata

    ds = utils.ensure_cf(ds)
    
    # Compute time stamps
    #ds = utils.shift_time(ds, ds.attrs["flowSampleDuration"] / 2)

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
                    iqmat[k][0:timelen, :], dims=("time", "bin_across")
                )
                if k in iqmat["Data_Units"]:
                    ds[k].attrs["units"] = iqmat["Data_Units"][k].replace("/s", " s-1")
            elif "_0_" in k or "_1_" in k:
                ds[k] = xr.DataArray(
                    iqmat[k][0:timelen, :], dims=("time", "bin_along")
                )
                if k in iqmat["Data_Units"]:
                    ds[k].attrs["units"] = iqmat["Data_Units"][k].replace("/s", " s-1")
            elif "FlowData_Vel" in k and "XYZ" not in k or "FlowData_SNR" in k:
                ds[k] = xr.DataArray(iqmat[k][0:timelen, :], dims=("time", "velbeam"))
                if k in iqmat["Data_Units"]:
                    ds[k].attrs["units"] = iqmat["Data_Units"][k].replace("/s", " s-1")
            elif "FlowData_VelXYZ" in k:
                ds["Vel_X_Center"] = xr.DataArray(iqmat[k][0:timelen, 0], dims=("time")
                )
                ds["Vel_Z_Center"] = xr.DataArray(iqmat[k][0:timelen, 1], dims=("time")
                )
                ds["Vel_X_Left"] = xr.DataArray(iqmat[k][0:timelen, 2], dims=("time")
                )
                ds["Vel_X_Right"] = xr.DataArray(iqmat[k][0:timelen, 3], dims=("time")
                )
                if k in iqmat["Data_Units"]:
                    xzvars = [
                        'Vel_X_Center',
                        'Vel_Z_Center',
                        'Vel_X_Left',
                        'Vel_X_Right'
                    ]
                    for var in xzvars:
                        ds[var].attrs["units"] = iqmat["Data_Units"][k].replace("/s", " s-1")
                    
            elif "FlowData_NoiseLevel" in k:
                ds[k] = xr.DataArray(iqmat[k][0:timelen, :], dims=("time", "beam"))
                if k in iqmat["Data_Units"]:
                    ds[k].attrs["units"] = iqmat["Data_Units"][k].replace("/s", " s-1")
                
        elif "FlowSubData" in k:
            if "CellSize" in k or "BlankingDistance" in k:
                ds[k] = (xr.DataArray(
                    iqmat[k][0:timelen], dims=("time")
                ) / 1000)
                if k in iqmat["Data_Units"]:
                    ds[k].attrs["units"] = iqmat["Data_Units"][k]
                
    ds["bin_along"] = np.arange(ds["Profile_0_Vel"].shape[1])
    ds["bin_across"] = np.arange(ds["Profile_2_Vel"].shape[1])

    ds = xr.Dataset(ds)

    for k in iqmat["System_IqSetup"]["basicSetup"]:
        if "spare" not in k:
            ds.attrs[k] = iqmat["System_IqSetup"]["basicSetup"][k]
    for k in iqmat["System_Id"]:
        ds.attrs[k] = iqmat["System_Id"][k]
    for k in iqmat["System_IqState"]:
        if "spare" not in k:
            ds.attrs[k.replace("[", "_").replace("]", "_")] = iqmat["System_IqState"][k]
            
    if ds.attrs["InstrumentType"] == "IQ":
        ds.attrs["AlongChannelBeamsAngle"] = 25
        ds.attrs["AcrossChannelBeamsAngle"] = 60
    else:
        print("Check and update beam angle for Sontek IQ instrument")

    return xr.decode_cf(ds)

def create_iqbindist(ds):
    """
    Generate bin distances from transducer along vertical profile for the along (beams 0,1) 
    and across (beams 2,3) bins. Cannot make bindist a dim because it changes by sample due 
    to bin size changing based on water depth. 
    """
    for bm in range(4):
        if bm < 2:
            bdname = "bin_along"
        else:
            bdname = "bin_across"

        r = range(len(ds[bdname]))
        cells = np.zeros(np.shape(ds["Profile_" + str(bm) + "_Vel"]))
        fsdbd = "FlowSubData_PrfHeader_" + str(bm) + "_BlankingDistance"
        fsdcs = "FlowSubData_PrfHeader_" + str(bm) + "_CellSize"
        for n in range(len(ds["time"])):
            #blanking distance + 1*bin_size = center of first bin
            #blanking distance + N*bin_size = center of N bin
            #due to 0 index, have to add 1 additional bin size to equation below
            #blanking distance + 1*binsize + (bin_along(or bin_across)*bin_size) = center of each bin
            cells[n, :] = ds[fsdbd][n].values + (r * ds[fsdcs][n].values) + (ds[fsdcs][n].values)
        ds["Profile_" + str(bm) + "_bindist"] = (xr.DataArray(cells, dims=("time", bdname)) )
        
        ds["Profile_" + str(bm) + "_bindist"].attrs.update(
            {
                "units": "m",
                "long_name": "bin (center) distance from transducer",
                "positive": "up",
                "note": "distance is along vertical profile from transducer",
            }
        )        

    return ds

def rename_vars(ds):
    # set up dict of instrument -> EPIC variable names
    
    newvars = {}
    
    for var in ds:
        if "Profile_0" in var:
            newvars[var] = var.replace("Profile_0_Amp","Profile_AGC1_1221").replace("Profile_0_Vel","Profile_vel1_1277").replace("Profile_0_BlankingDistance","Profile_blanking_distance1").replace("Profile_0_CellSize","Profile_bin_size1").replace("Profile_0_bindist","Profile_bindist1").replace("Profile_0_z","Profile_z1").replace("Profile_0_bindepth","Profile_bindepth1")
        elif "Profile_1" in var:
            newvars[var] = var.replace("Profile_1_Amp","Profile_AGC2_1222").replace("Profile_1_Vel","Profile_vel2_1278").replace("Profile_1_BlankingDistance","Profile_blanking_distance2").replace("Profile_1_CellSize","Profile_bin_size2").replace("Profile_1_bindist","Profile_bindist2").replace("Profile_1_z","Profile_z2").replace("Profile_1_bindepth","Profile_bindepth2")
        elif "Profile_2" in var:
            newvars[var] = var.replace("Profile_2_Amp","Profile_AGC3_1223").replace("Profile_2_Vel","Profile_vel3_1279").replace("Profile_2_BlankingDistance","Profile_blanking_distance3").replace("Profile_2_CellSize","Profile_bin_size3").replace("Profile_2_bindist","Profile_bindist3").replace("Profile_2_z","Profile_z3").replace("Profile_2_bindepth","Profile_bindepth3")
        elif "Profile_3" in var:
            newvars[var] = var.replace("Profile_3_Amp","Profile_AGC4_1224").replace("Profile_3_Vel","Profile_vel4_1280").replace("Profile_3_BlankingDistance","Profile_blanking_distance4").replace("Profile_3_CellSize","Profile_bin_size4").replace("Profile_3_bindist","Profile_bindist4").replace("Profile_3_z","Profile_z4").replace("Profile_3_bindepth","Profile_bindepth4")
            
    
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
    for k in varnames:
        if k in ds:
            newvars[k] = varnames[k]
        
    return ds.rename(newvars)


def update_prefixes(ds):
    newvars = {}
    for k in ds:
        if "FlowData" in k:
            newvars[k] = k.replace("FlowData_", "")
        
        elif "FlowSubData_PrfHeader" in k:
                newvars[k] = k.replace("FlowSubData_PrfHeader_", "Profile_")            

    return ds.rename(newvars)

def vel_to_ms(ds):
    """
    Convert velocity data from mm/s to m/s
    """

    for var in ds:
        if "Vel" in var:
            ds[var] = ds[var] / 1000
            ds[var].attrs['units'] = "m s-1"

    return ds

def create_iqbindepth(ds):
    """
    Generate bin depths reltive to pressure. 
    """
    for bm in range(4):
        ds['Profile_%d_bindepth' % bm] = ds['AdjustedPressure'] - ds['Profile_%d_bindist' % bm]
        ds["Profile_%d_bindepth" % bm].attrs.update(
            {
                "units": "m",
                "long_name": "bin(center) depth relative to sea surface",
                "positive": "down",
                "note": "Distance is along vertical profile from transducer",
            }
        )
        
    return ds

def create_iqz(ds):
    """
    Generate bin heights relative to geopotential datum
    """
    for bm in range(4):
        if "height_above_geopotential_datum" in ds.attrs:
            if ds.attrs["orientation"].upper() == "DOWN":
                if bm < 2:
                    ds['Profile_%d_z' % bm] = xr.DataArray(
                        ds.attrs["height_above_geopotential_datum"]
                        + ds.attrs["initial_instrument_height"]
                        - ds["Profile_%d_bindist" % bm].values,
                        dims=("time", "bin_along"),
                    )
                else:
                    ds['Profile_%d_z' % bm] = xr.DataArray(
                        ds.attrs["height_above_geopotential_datum"]
                        + ds.attrs["initial_instrument_height"]
                        - ds["Profile_%d_bindist" % bm].values,
                        dims=("time", "bin_across"),
                    )                    
            elif ds.attrs["orientation"].upper() == "UP":
                if bm < 2:
                    ds['Profile_%d_z' % bm] = xr.DataArray(
                        ds.attrs["height_above_geopotential_datum"]
                        + ds.attrs["initial_instrument_height"]
                        + ds["Profile_%d_bindist" % bm].values,
                        dims=("time", "bin_along"),
                    )
                else:
                    ds['Profile_%d_z' % bm] = xr.DataArray(
                        ds.attrs["height_above_geopotential_datum"]
                        + ds.attrs["initial_instrument_height"]
                        + ds["Profile_%d_bindist" % bm].values,
                        dims=("time", "bin_across"),
                    )                 
                        
            else:
                print("Could not create z for bins, specifiy orientation")
                
        else:
            print("Could not create z for bins, specify height_above_geopotential_datum")
            
    return (ds)

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
    Load a raw .cdf file and generate a processed .nc file
    """

    # Load raw .cdf data
    ds = xr.open_dataset(cdf_filename)

    ds = update_prefixes(ds)

    # Clip data to in/out water times or via good_ens
    ds = utils.clip_ds(ds)

    ds = vel_to_ms(ds)
    
    ds = create_iqbindepth(ds)
    
    ds = create_iqz(ds)

    ds = clean_iq(ds)
    
    ds = trim_iqvel(ds)
    
    ds = fill_snr(ds)
    
    ds = fill_vbper(ds)

    # assign min/max:
    ds = utils.add_min_max(ds)

    ds = utils.add_start_stop_time(ds)

    ds = utils.add_delta_t(ds)

    # add lat/lon coordinates
    ds = utils.ds_add_lat_lon(ds)

    ds = rename_vars(ds)
           
    # should function this
    for var in ds.data_vars:
        ds = qaqc.trim_min(ds, var)
        ds = qaqc.trim_max(ds, var)
        ds = qaqc.trim_min_diff(ds, var)
        ds = qaqc.trim_max_diff(ds, var)
        ds = qaqc.trim_med_diff(ds, var)
        ds = qaqc.trim_med_diff_pct(ds, var)
        ds = qaqc.trim_bad_ens(ds, var)
        ds = qaqc.trim_maxabs_diff_2d(ds, var)
        ds = qaqc.trim_fliers(ds, var)

    # after check for masking vars by other vars
    for var in ds.data_vars:
        ds = qaqc.trim_mask(ds, var)
        
    ds = utils.create_z(ds) #added 7/31/2023

    ds = ds_add_attrs(ds)

    #ds = utils.no_p_create_depth(ds) #commented out 7/31/23

    ds = ds.drop(["SampleNumber", "SampleTime","Volume_Total","Volume_Positive","Volume_Negative"])

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
    
    newvars = {}
    for k in dsprof:
        newvars[k] = k.replace("Profile_", "")
        
    dsprof = dsprof.rename(newvars)    

    # Write to .nc file
    print("Writing cleaned/trimmed data to .nc file")

    nc_filename = dsflow.attrs["filename"] + "flow-a.nc"
    dsflow.to_netcdf(nc_filename, unlimited_dims=["time"])
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])
    print("Done writing netCDF file", nc_filename)

    nc_filename = dsprof.attrs["filename"] + "prof-a.nc"
    dsprof.to_netcdf(nc_filename)
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

    ds["bin_across"].attrs.update(
        {"long_name": "bin number for across channel/skew beams (beams 3 & 4)", "units": "1", "bin_size": "size of bin varies by sample. See variables bin_size3 and bin_size4", "bin_count": "%d" % len(ds.bin_across), "blanking_distance": "blanking distance varies by sample. See profile variables blanking_distance3 and blanking_distance4", "note": "bin number is along profile from corresponding transducer"}
    )
    
    ds["bin_along"].attrs.update(
        {"long_name": "bin number for along channel beams (beams 1 & 2)", "units": "1", "bin_size": "size of bin varies by sample. See variables bin_size1 and bin_size2", "bin_count": "%d" % len(ds.bin_across), "blanking_distance": "blanking distance varies by sample. See profile variables blanking_distance1 and blanking_distance2", "note": "bin number is along profile from corresponding transducer"}
    ) 
    
    ds["velbeam"].attrs.update(
        {"long_name": "velocity beam number", "units": "1", "note":"does not include vertical beam (vb, beam 5)"}
    )
 
    ds["beam"].attrs.update(
        {"long_name": "beam number", "units": "1", "note":"includes vertical beam (vb, beam 5)"}
    )         
         
    ds["D_3"].attrs.update(
        {"long_name": "acoustically measured depth", "note": "relative to vertical transducer/beam 5, the top of the instrument"}
    )

    # descriptions from Sontek-IQ Series User's Manual available at
    # http://info.xylem.com/sontek-iq-manual.html
    ds["Stage"].attrs["long_name"] = "Stage (water depth of the user-defined channel)"
    ds["Area"].attrs["long_name"] = "Cross-sectional area of user-defined channel"
    ds["Flow"].attrs["long_name"] = "Flow rate (using defined channel geometry)"

    ds["Vel_Mean"].attrs.update({"long_name" : "Mean velocity", "positive_dir" : "%s" % ds.attrs["beam_2_positive_direction"]})
    ds["Volume_Total"].attrs[
        "long_name"
    ] = "Total water volume (based on all measured flow)"
    ds["Volume_Positive"].attrs[
        "long_name"
    ] = "Total volume of water in the positive downstream direction"
    ds["Volume_Negative"].attrs[
        "long_name"
    ] = "Total volume of water in the negative upstream direction"
    ds["Vel"].attrs.update({"long_name" : "Beam velocity", "beam_2_positive_dir" : "%s" % ds.attrs['beam_2_positive_direction'], "beams_1,3,4_positive_dir" : "%s" % ds.attrs['beam_1_positive_direction']})
    ds["Vel_X_Center"].attrs.update({"long_name" : "X velocity from center beams (beams 1 & 2)", "positive_dir" : "%s" % ds.attrs['beam_2_positive_direction']})
    ds["Vel_Z_Center"].attrs.update({"long_name" : "Z velocity from center beams (beams 1 & 2)", "positive_dir" : "%s" % ds.attrs['orientation']})
    ds["Vel_X_Left"].attrs.update({"long_name" : "X velocity left beam (beam 3)", "positive_dir" : "%s" % ds.attrs['beam_2_positive_direction']})
    ds["Vel_X_Right"].attrs.update({"long_name" : "X velocity right beam (beam 4)", "positive_dir" : "%s" % ds.attrs['beam_2_positive_direction']})
    ds["VelStd"].attrs["long_name"] = "Velocity standard deviation"
    ds["SNR"].attrs["long_name"] = "Signal-to-noise ratio"
    ds["NoiseLevel"].attrs["long_name"] = "Acoustic noise level"
    ds["Range"].attrs["long_name"] = "Acoustically measured distance to water surface"
    ds["T_28"].attrs.update({"long_name" : "Temperature", "epic_code":"28", "units":"degree_C","standard_name":"sea_water_temperature"})
    ds["P_1"].attrs.update({"long_name" : "Uncorrected pressure", "epic_code":"1", "standard_name":"sea_water_pressure"})
    ds["PressOffsetAdjust"].attrs.update({"long_name":"Atmospheric pressure adjustment", "note":"see SonTek-IQ User's Manual for details"})
    ds["P_1ac"].attrs.update({"long_name" : "Corrected pressure", "standard_name" : "sea_water_pressure_due_to_sea_water", "note" : "Measurement with atmospheric pressure removed (see SonTek-IQ User's Manual for details)"})
    ds["Bat_106"].attrs.update({"long_name" : "Battery voltage", "epic_code":"106"})
    ds["Ptch_1216"].attrs["long_name"] = "Pitch angle in degrees"
    
    # to be UDUNITS compatible
    if ds["Ptch_1216"].attrs["units"] == "deg":
        ds["Ptch_1216"].attrs["units"] = "degrees"
    if ds["Roll_1217"].attrs["units"] == "deg":
        ds["Roll_1217"].attrs["units"] = "degrees"
    ds["Roll_1217"].attrs.update({"long_name":"Instrument Roll","epic_code":"1217","standard_name":"platform_roll"})
    ds["Ptch_1216"].attrs.update({"long_name":"Instrument Pitch","epic_code":"1216","standard_name":"platform_pitch"})
    ds["VbPercentGood"].attrs["long_name"] = "Vertical beam percent good"
    ds["HorizontalSkew"].attrs["long_name"] = "Horizontal skew"
    ds["SystemInWater"].attrs.update({"long_name" : "Percentage of sample during which instrument was submerged", "note":"100% means it was submerged for entire sample"})

    # Profile Variables
    for n in range(4):
        ds["Profile_AGC%d_122%d" % (n + 1, n +1)].attrs.update({"units":"counts","long_name" : "Echo Intensity (AGC) beam %d" % (n + 1)})
        ds["Profile_vel%d_%dStd" % (n + 1, n + 1277)].attrs["long_name"] = (
            "beam %d velocity profile standard deviation" % (n + 1)
        )
        ds["Profile_vel%d_%d" % (n + 1, n + 1277)].attrs.update({"long_name" : "beam %d current velocity" % (n + 1), "positive_dir" : "%s" % ds.attrs["beam_" + str(n+1) + "_positive_direction"]})
        ds["Profile_blanking_distance%d" % (n + 1)].attrs.update({"long_name":"beam %d blanking distance" % (n + 1), "units": "m"})
        ds["Profile_bin_size%d" % (n + 1)].attrs.update({"long_name":"beam %d bin size" % (n + 1), "units": "m"}) 
        ds["Profile_z%d" % (n + 1)].attrs.update({"standard_name": "height", "long_name":"beam %d bin height relative to %s" % ((n + 1), ds.attrs["geopotential_datum_name"]), "units": "m", "positive":"%s" % ds.attrs["orientation"], "axis" : "Z"})
         
         
    return ds

def trim_iqvel(ds):
    """
    Trim velocity data depending on specified method
    """
        
    if (
        "trim_method" in ds.attrs
        and ds.attrs["trim_method"].lower() != "none"
        and ds.attrs["trim_method"] is not None
    ):        
        if "AdjustedPressure" in ds:
            P = ds["AdjustedPressure"]
            Ptxt = "atmospherically corrected"
        elif "P_1ac" in ds:
            P = ds["P_1ac"]
            Ptxt = "atmospherically corrected"
        elif "Pressure" in ds:
            # FIXME incorporate press_ ac below
            P = ds["Pressure"]
            Ptxt = "NON-atmospherically corrected"
            
        for bm in range(4): #beams are different angles
            if bm < 2:
                bmangle = ds.attrs["AlongChannelBeamsAngle"]
            else:
                bmangle = ds.attrs["AcrossChannelBeamsAngle"]
                
            if ds.attrs["trim_method"].lower() == "water level":
                ds["Profile_" + str(bm) + "_Vel"] = ds["Profile_" + str(bm) + "_Vel"].where(ds["Profile_" + str(bm) + "_bindist"] < P)
            
                histtext = "Trimmed velocity data using {} pressure (water level).".format(
                    Ptxt
                )

                
            elif ds.attrs["trim_method"].lower() == "water level sl":
                ds["Profile_" + str(bm) + "_Vel"] = ds["Profile_" + str(bm) + "_Vel"].where(
                ds["Profile_" + str(bm) + "_bindist"] < P * np.cos(np.deg2rad(bmangle))
                )

                histtext = "Trimmed velocity data using {} pressure (water level) and sidelobes.".format(
                    Ptxt
                )

        ds = utils.insert_history(ds, histtext)
            
    else:
        print("Did not trim velocity data")
        
    return(ds)

def fill_snr(ds):
    """
    Fill velocity data with corresponding beam snr value threshold
    """
    if "snr_threshold" in ds.attrs:
        
        Ptxt = str(ds.attrs['snr_threshold'])
    
        for var in ds:
            if "Vel" and "Profile" in var:
                for bm in range(4):
                    var = "Profile_" + str(bm) + "_Vel"
                    ds[var] = ds[var].where(ds.SNR[:,bm] > ds.attrs['snr_threshold'])
    
            else:
                ds["Vel"] = ds["Vel"].where(ds.SNR > ds.attrs['snr_threshold'])
                ds["Vel_X_Center"] = ds["Vel_X_Center"].where((ds.SNR[:,0] > ds.attrs['snr_threshold']) & (ds.SNR[:,1] > ds.attrs['snr_threshold']))
                ds["Vel_Z_Center"] = ds["Vel_Z_Center"].where((ds.SNR[:,0] > ds.attrs['snr_threshold']) & (ds.SNR[:,1] > ds.attrs['snr_threshold']))
                ds["Vel_X_Left"] = ds["Vel_X_Left"].where(ds.SNR[:,2] > ds.attrs['snr_threshold'])
                ds["Vel_X_Right"] = ds["Vel_X_Right"].where(ds.SNR[:,3] > ds.attrs['snr_threshold'])
                ds["Vel_Mean"] = ds["Vel_Mean"].where((ds.SNR[:,0] > ds.attrs['snr_threshold']) & (ds.SNR[:,1] > ds.attrs['snr_threshold']) & (ds.SNR[:,2] > ds.attrs['snr_threshold']) & (ds.SNR[:,3] > ds.attrs['snr_threshold']))
            
            histtext = "Filled velocity data using snr threshold of {} for corresponding beams.".format(
                Ptxt
            )

        ds = utils.insert_history(ds, histtext)
    else:
            print('Did not fill velocity data using snr threshold')
    return(ds)

def fill_vbper(ds):
    """
    Fill adjusted pressure, stage, area, range, and depth data with corresponding vertical beam percent good threshold
    """

    if "vbper_threshold" in ds.attrs:
        
        Ptxt = str(ds.attrs['vbper_threshold'])
        
        histtext = "Filling P1ac, stage, area, range, D_3, and profile velocity data using vertical beam percent good threshold threshold of {}.".format(Ptxt)
        
        for var in ds:
            if "Vel" and "Profile" in var:
                for bm in range(4):
                    var = "Profile_" + str(bm) + "_Vel"
                    ds[var] = ds[var].where(ds.VbPercentGood > ds.attrs['vbper_threshold'])
                    
            else:
                
                ds["AdjustedPressure"] = ds["AdjustedPressure"].where(ds.VbPercentGood > ds.attrs['vbper_threshold'])
                ds["Depth"] = ds["Depth"].where(ds.VbPercentGood > ds.attrs['vbper_threshold'])
                ds["Stage"] = ds["Stage"].where(ds.VbPercentGood > ds.attrs['vbper_threshold'])
                ds["Area"] = ds["Area"].where(ds.VbPercentGood > ds.attrs['vbper_threshold'])
                ds["Range"] = ds["Range"].where(ds.VbPercentGood > ds.attrs['vbper_threshold'])
                ds["Vel_Mean"] = ds["Vel_Mean"].where(ds.VbPercentGood > ds.attrs['vbper_threshold'])
            
        ds = utils.insert_history(ds, histtext)
        
    else:
        
        print('Did not fill pressure, stage, area, range, and depth data data using snr threshold')
            
    return(ds)        