import pkgutil
import time

import xarray as xr

#import dolfyn
import xmltodict

from ..aqd import aqdutils
from ..core import utils


def cdf_to_nc(cdf_filename, atmpres=False):
    """
    Load a "raw" .cdf file and generate a processed .nc file
    """
    print(f"Loading {cdf_filename[0]}")
    start_time = time.time()
    ds=xr.open_dataset(cdf_filename,chunks="auto")
    end_time = time.time()
    print(f"Finished loading {cdf_filename[0]} in {end_time-start_time:.1f} seconds")


    """
    # TODO: Add atmospheric pressure offset
    print(f"Loading {cdf_filename[0]}")
    start_time = time.time()
    ds = dolfyn.load(cdf_filename[0])
    end_time = time.time()
    print(f"Finished loading {cdf_filename[0]} in {end_time-start_time:.1f} seconds")
    """
    ds = utils.create_nominal_instrument_depth(ds)

    # create Z depending on orientation
    #ds, T, T_orig = aqdutils.set_orientation(ds, ds["Burst_Beam2xyz"].values)
    ds=utils.create_z(ds)

    """# Transform ("rotate") into ENU coordinates
    dolfyn.rotate2(ds, "earth")"""

    # Clip data to in/out water times or via good_ens
    # Need to clip data after coord transform when using dolfyn
    ds = utils.clip_ds(ds)

    """# Create separate vel variables first
    ds["U"] = ds["vel"].sel(dir="E")
    ds["V"] = ds["vel"].sel(dir="N")
    ds["W1"] = ds["vel"].sel(dir="U1")
    ds["W2"] = ds["vel"].sel(dir="U2")"""
    
    ds["U"] = ds["VelEast"]
    ds["V"] = ds["VelNorth"]
    ds["W1"] = ds["VelUp1"]
    ds["W2"] = ds["VelUp2"]

    ds = aqdutils.magvar_correct(ds)
    
    ds = aqdutils.trim_vel(ds,data_vars=["U","V","W1","W2"])
    ds = aqdutils.make_bin_depth(ds)

    # Rename DataArrays for EPIC compliance
    ds = aqdutils.ds_rename(ds) #for common variables
    ds = ds_rename_sig(ds) #for signature vars not in aqds or vecs
    ds=ds_drop(ds)
    # swap_dims from bindist to depth
    ds = ds_swap_dims(ds)

    # Add EPIC and CMG attributes
    ds = aqdutils.ds_add_attrs(ds, inst_type="SIG")

     # Add min/max values
    ds = utils.add_min_max(ds)

    # Add DELTA_T for EPIC compliance
    #ds = utils.add_delta_t(ds)

    # Add start_time and stop_time attrs
    ds = utils.add_start_stop_time(ds)

    # Add history showing file used
    ds = utils.add_history(ds)

    """ds = clean_dolfyn_standard_names(ds)"""

    ds = utils.add_standard_names(ds)

    """ds = fix_dolfyn_encoding(ds)"""

    """# split up into multiple files and write them separately
    ds_b5 = ds.copy()
    ds_echo = ds.copy()"""

    todrop = []
    """for v in ds.data_vars:
        if "_b5" in v or "_echo" in v:
            todrop.append(v)"""

    ds = ds.drop_vars(todrop)
    ds = drop_unused_dims(ds)

    """todrop = []
    for v in ds_b5.data_vars:
        if "_b5" not in v:
            todrop.append(v)

    ds_b5 = ds_b5.drop_vars(todrop)
    ds_b5 = drop_unused_dims(ds_b5)

    todrop = []
    for v in ds_echo.data_vars:
        if "_echo" not in v:
            todrop.append(v)

    ds_echo = ds_echo.drop_vars(todrop)
    ds_echo = drop_unused_dims(ds_echo)"""

    """if "prefix" in ds.attrs:
        nc_filename = ds.attrs["prefix"] + ds.attrs["filename"] + "-a.nc"
    else:
        nc_filename = ds.attrs["filename"] + "-a.nc"

    for datatype, dsout in zip(["", "_b5", "_echo"], [ds, ds_b5, ds_echo]):
        if datatype == "":
            nc_out = nc_filename
        elif datatype == "_b5":
            nc_out = nc_filename[:-5] + "_b5-a.nc"
        elif datatype == "_echo":
            nc_out = nc_filename[:-5] + "_echo-a.nc"

        dolfyn.save(dsout, nc_out)
        utils.check_compliance(nc_out, conventions=ds.attrs["Conventions"])

        print("Done writing netCDF file", nc_filename)"""
    #write out nc file by data_type
    if "prefix" in ds.attrs:
        nc_filename = ds.attrs["prefix"] + ds.attrs["filename"]
    else:
        nc_filename = ds.attrs["filename"]

    if ds.attrs["data_type"] == "Burst":
        if "sample" not in ds.dims:
            ds["time"].encoding["dtype"] = "double"
        else:
            ds["time"].encoding["dtype"] = "i4"
        
        nc_out= nc_filename + "b-cal.nc"
        print("writing Burst (b) data to netCDF nc file")
        ds.to_netcdf(nc_out)
        print(f"Finished writing data to {nc_out}")

    elif ds.attrs["data_type"] == "IBurst":
        if "sample" not in ds.dims:
            ds["time"].encoding["dtype"] = "double"
        else:
            ds["time"].encoding["dtype"] = "i4"
        nc_out= nc_filename + "b5-cal.nc"
        print("writing IBurst (b5) data to netCDF nc file")
        ds.to_netcdf(nc_out)
        print(f"Finished writing data to {nc_out}")
        

    utils.check_compliance(nc_out, conventions=ds.attrs["Conventions"])

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


def clean_dolfyn_standard_names(ds):
    """remove non-compliant standard_names set by dolfyn"""
    data = pkgutil.get_data(__name__, "../data/cf-standard-name-table.xml")
    doc = xmltodict.parse(data)

    entries = doc["standard_name_table"]["entry"]
    allnames = [x["@id"] for x in entries]

    for v in ds.data_vars:
        if "standard_name" in ds[v].attrs:
            if ds[v].attrs["standard_name"] not in allnames:
                del ds[v].attrs["standard_name"]

    return ds


def fix_dolfyn_encoding(ds):
    """ensure we don't set dtypes of int64 for CF compliance"""
    for var in ds.dims:
        if ds[var].dtype == "int64":
            if ds[var].max() > 2**31 - 1 or ds[var].min() < -(2**31):
                print(
                    f"warning {var} may be too big to fit in int32: min {ds[var].min().values}, max {ds[var].max().values}"
                )
            ds[var].encoding["dtype"] = "int32"

    return ds

def ds_drop(ds):
    """
    Drop old DataArrays from Dataset that won't make it into the final .nc file
    """

    todrop = [
        "ExtStatus",
        "NBeams",
        "NCells",
        "PressureSensorTemperature",
        "RTCTemperature",
        "MagnetometerTemperature",
        "VEL1",
        "VEL2",
        "VEL3",
        "AMP1",
        "AMP2",
        "AMP3",
        "VelEast",
        "VelNorth",
        "VelUp1",
        "VelUp2",
        "VelX",
        "VelY",
        "VelZ1",
        "VelZ2",
        "Beam2xyz"
    ]

    if ("AnalogInput1" in ds.attrs) and (ds.attrs["AnalogInput1"].lower() == "true"):
        todrop.remove("AnalogInput1")

    if ("AnalogInput2" in ds.attrs) and (ds.attrs["AnalogInput2"].lower() == "true"):
        todrop.remove("AnalogInput2")

    return ds.drop([t for t in todrop if t in ds.variables])

def ds_rename_sig(ds, waves=False):
    """
    Rename DataArrays within Dataset for compliance
    """    
    varnames= (
        {
            "EnsembleCount":"sample",
            "AmbiguityVel":"AmbVel",
            "U": "u_1205",
            "V": "v_1206",
            "W1": "w_1204",
            "W2": "w2_1204",
            "VelSpeed":"CS_300",
            "VelDirection":"CD_310",
            "VelBeam1":"vel1_1277",
            "VelBeam2":"vel1_1278",
            "VelBeam3":"vel1_1279",
            "VelBeam4":"vel1_1280",
            "AmpBeam1":"Sv1_1233",
            "AmpBeam2":"Sv2_1234",
            "AmpBeam3":"Sv3_1235",
            "AmpBeam4":"Sv1_1236",
            "CorBeam1":"cor1_1285",
            "CorBeam2":"cor2_1286",
            "CorBeam3":"cor3_1287",
            "CorBeam4":"cor4_1288",
            "AHRSRotationMatrix":"orientmat",
        }
    )

    for v in varnames:
        if v in ds:
            ds = ds.rename({v: varnames[v]})

    for v in [
        "avgamp1",
        "avgamp2",
        "avgamp3",
        "U",
        "V",
        "W",
        "Depth",
        "water_depth",
        "cellpos",
        "vel",  # Signature velocity
    ]:
        if v in ds:
            ds = ds.drop_vars(v)

    return ds

def ds_swap_dims(ds):
    # need to preserve z attrs because swap_dims will remove them
    attrsbak = ds["z"].attrs
    for v in ds.data_vars:
        if "bindist" in ds[v].coords:
            ds[v] = ds[v].swap_dims({"bindist": "z"})

    ds["z"].attrs = attrsbak

    return ds