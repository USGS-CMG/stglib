import xarray as xr
import numpy as np

from ..core import utils

from ..aqd import aqdutils


def cdf_to_nc(cdf_filename, atmpres=None):
    """
    Load a "raw" .cdf file and generate a processed .nc file
    """

    # Load raw .cdf data
    ds = xr.load_dataset(cdf_filename)

    # Clip data to in/out water times or via good_ens
    ds = utils.clip_ds(ds)

    ds = utils.create_nominal_instrument_depth(ds)

    ds = change_units(ds)

    ds = create_depth(ds)

    if atmpres is not None:
        print("Atmospherically correcting data")

        met = xr.load_dataset(atmpres)
        # need to save attrs before the subtraction, otherwise they are lost
        # ds['P_1ac'] = ds['P_1'].copy(deep=True)
        attrs = ds["P_1"].attrs
        ds["P_1ac"] = ds["P_1"] - met["atmpres"] - met["atmpres"].offset
        print("Correcting using offset of %f" % met["atmpres"].offset)
        ds["P_1ac"].attrs = attrs

    if ds.attrs["RDI_Coord_Transform"] == "BEAM":
        rotmat = calc_beam_rotmatrix(
            ds.attrs["RDI_Beam_Angle"], ds.attrs["RDI_Beam_Pattern"].lower() == "convex"
        )

        inst = np.full((ds.vel1.shape[0], ds.vel1.shape[1], 4), np.nan)
        earth = np.full((ds.vel1.shape[0], ds.vel1.shape[1], 4), np.nan)
        print("Converting from beam to instrument coordinates")
        for n in range(len(ds.vel1[:, 0])):
            vels = np.vstack(
                [ds.vel1[n, :], ds.vel2[n, :], ds.vel3[n, :], ds.vel4[n, :]]
            ).T
            tmp = beam2inst(vels, rotmat=rotmat)
            inst[n, :, :] = beam2inst(vels, rotmat=rotmat)
            # print(inst[n,:,:].shape)
        print("Converting from instrument to earth coordinates")
        for n in range(len(ds.vel1[:, 0])):
            earth[n, :, :] = inst2earth(
                inst[n, :, :],
                ds["Hdg"][n].values,
                ds["Ptch"][n].values,
                ds["Roll"][n].values,
            )

    for n, var in enumerate(["u_1205", "v_1206", "w_1204", "Werr_1201"]):
        ds[var] = xr.DataArray(earth[:, :, n], dims=("time", "bindist"))
        ds[var].attrs["units"] = ds["vel1"].attrs["units"]

    ds["u_1205"].attrs["standard_name"] = "eastward_sea_water_velocity"
    ds["v_1206"].attrs["standard_name"] = "northward_sea_water_velocity"
    ds["w_1204"].attrs["standard_name"] = "upward_sea_water_velocity"
    ds["Werr_1201"].attrs[
        "standard_name"
    ] = "indicative_error_from_multibeam_acoustic_doppler_velocity_profiler_in_sea_water"

    ds = aqdutils.magvar_correct(ds)

    for todrop in ["Pressure+", "Pressure-", "EWD1", "EWD2", "EWD3", "EWD4"]:
        ds = ds.drop(todrop)

    ds["time"].attrs["standard_name"] = "time"
    ds["time"].encoding["dtype"] = "int32"

    # get dtypes right for cf compliance
    for var in ds:
        if ds[var].dtype == "uint16":
            ds[var].encoding["dtype"] = "int32"
        elif ds[var].dtype == "int64":
            if ds[var].max() > 2 ** 31 - 1 or ds[var].min() < -(2 ** 31):
                print(
                    f"warning {var} may be too big to fit in int32: min {ds[var].min().values}, max {ds[var].max().values}"
                )
            ds[var].encoding["dtype"] = "int32"

    # clean up attrs for cf compliance letters digits and underscores only.
    for attr in ds.attrs:
        if "." in attr:
            ds.attrs[attr.replace(".", "_")] = ds.attrs.pop(attr)
        if "-" in attr:
            ds.attrs[attr.replace("-", "_")] = ds.attrs.pop(attr)
        if " " in attr:
            ds.attrs[attr.replace(" ", "_")] = ds.attrs.pop(attr)
        if "+" in attr:
            ds.attrs[attr.replace(" ", "_")] = ds.attrs.pop(attr)

    ds.to_netcdf(ds.attrs["filename"] + "-clean.nc")
    return ds

    # Create depth variable depending on orientation
    VEL, T, T_orig = aqdutils.set_orientation(VEL, VEL["TransMatrix"].values)

    # Transform coordinates from, most likely, BEAM to ENU
    u, v, w = aqdutils.coord_transform(
        VEL["VEL1"].values,
        VEL["VEL2"].values,
        VEL["VEL3"].values,
        VEL["Heading"].values,
        VEL["Pitch"].values,
        VEL["Roll"].values,
        T,
        T_orig,
        VEL.attrs["AQDCoordinateSystem"],
    )

    VEL["U"] = xr.DataArray(u, dims=("time", "bindist"))
    VEL["V"] = xr.DataArray(v, dims=("time", "bindist"))
    VEL["W"] = xr.DataArray(w, dims=("time", "bindist"))

    VEL = aqdutils.magvar_correct(VEL)

    VEL["AGC"] = (VEL["AMP1"] + VEL["AMP2"] + VEL["AMP3"]) / 3

    VEL = aqdutils.trim_vel(VEL)

    VEL = aqdutils.make_bin_depth(VEL)

    # Reshape and associate dimensions with lat/lon
    for var in [
        "U",
        "V",
        "W",
        "AGC",
        "Pressure",
        "Temperature",
        "Heading",
        "Pitch",
        "Roll",
        "bin_depth",
        "Pressure_ac",
    ]:
        if var in VEL:
            VEL = utils.add_lat_lon(VEL, var)

    # swap_dims from bindist to depth
    VEL = ds_swap_dims(VEL)

    # Rename DataArrays for EPIC compliance
    VEL = aqdutils.ds_rename(VEL)

    # Drop non-EPIC variables
    VEL = ds_drop(VEL)

    # Add EPIC and CMG attributes
    VEL = aqdutils.ds_add_attrs(VEL)

    # Add min/max values
    VEL = utils.add_min_max(VEL)

    # Add DELTA_T for EPIC compliance
    VEL = aqdutils.add_delta_t(VEL)

    # Add start_time and stop_time attrs
    VEL = utils.add_start_stop_time(VEL)

    # Add history showing file used
    VEL = utils.add_history(VEL)

    if utils.is_cf(VEL):
        VEL = utils.add_standard_names(VEL)

    for var in VEL.variables:
        if (var not in VEL.coords) and ("time" not in var):
            # cast as float32
            VEL = utils.set_var_dtype(VEL, var)

    if "prefix" in VEL.attrs:
        nc_filename = VEL.attrs["prefix"] + VEL.attrs["filename"] + "-a.nc"
    else:
        nc_filename = VEL.attrs["filename"] + "-a.nc"

    if utils.is_cf(VEL):
        VEL.to_netcdf(nc_filename, encoding={"time": {"dtype": "i4"}})
    else:
        VEL.to_netcdf(nc_filename, unlimited_dims=["time"])

    print("Done writing netCDF file", nc_filename)

    return VEL


# TODO: add analog input variables (OBS, NTU, etc)


def create_depth(ds):
    """
    Create depth variable depending on instrument orientation
    """

    if "Pressure_ac" in ds:
        presvar = "Pressure_ac"
    elif "Pressure" in ds:
        presvar = "Pressure"
    else:
        presvar = None

    if "NAVD88_ref" in ds.attrs:
        Wdepth = -ds.attrs["NAVD88_ref"] - ds.attrs["transducer_offset_from_bottom"]
    elif presvar is not None:
        Wdepth = np.nanmean(ds[presvar])
    else:
        Wdepth = ds.attrs["nominal_instrument_depth"]

    if ds.attrs["orientation"] == "UP":
        print("User instructed that instrument was pointing UP")
        # try a simpler approach
        ds["depth"] = xr.DataArray(Wdepth - ds["bindist"], dims="bindist")
    # elif ds.attrs["orientation"] == "DOWN":
    #     print("User instructed that instrument was pointing DOWN")
    #     T[1, :] = -T[1, :]
    #     T[2, :] = -T[2, :]
    #     # try a simpler approach
    #     VEL["depth"] = xr.DataArray(Wdepth + VEL["bindist"], dims="bindist")

    return ds


def ds_swap_dims(ds):
    ds = ds.swap_dims({"bindist": "depth"})
    # need to swap dims and then reassign bindist to be a normal variable
    # (no longer a coordinate)
    valbak = ds["bindist"].values
    ds = ds.drop("bindist")
    ds["bindist"] = xr.DataArray(valbak, dims="depth")

    return ds


def ds_drop(ds):
    """
    Drop old DataArrays from Dataset that won't make it into the final .nc file
    """

    todrop = [
        "VEL1",
        "VEL2",
        "VEL3",
        "AMP1",
        "AMP2",
        "AMP3",
        "TransMatrix",
        "AnalogInput1",
        "AnalogInput2",
        "jd",
        "Depth",
    ]

    return ds.drop([t for t in todrop if t in ds.variables])


def change_units(ds):

    # pressure needs to be in db or m
    if "Pressure" in ds:
        pconvconst = 1  # when in doubt, do nothing
        if "deca-pascals" in ds["Pressure"].attrs["units"]:
            pconvconst = 1000  # decapascals to dbar = /1000
            print("Pressure in deca-pascals will be converted to db")
            ds["Pressure"].attrs["units"] = "dbar"
        attrsbak = ds["Pressure"].attrs
        ds["Pressure"] = ds["Pressure"] / pconvconst
        for k in attrsbak:
            ds["Pressure"].attrs[k] = attrsbak[k]

    # check units of current velocity and convert to cm/s
    vunits = ds["vel1"].attrs["units"]
    vconvconst = 1  # when in doubt, do nothing
    if (vunits == "mm s-1") | (vunits == "mm/s"):
        vconvconst = 0.1  # mm/s to cm/s
    elif (vunits == "m s-1") | (vunits == "m/s"):
        vconvconst = 100  # m/s to cm/s
    for var in ["vel1", "vel2", "vel3", "vel4"]:
        ds[var].attrs["units"] = "cm s-1"
        attrsbak = ds[var].attrs
        ds[var] = ds[var] * vconvconst
        for k in attrsbak:
            ds[var].attrs[k] = attrsbak[k]

    for var in ds:
        if ds[var].attrs["units"] == "hundredths of degrees":
            ds[var].attrs["units"] = "degree"
            attrsbak = ds[var].attrs
            ds[var] = ds[var] / 100
            for k in attrsbak:
                ds[var].attrs[k] = attrsbak[k]

    for var in ds:
        if ds[var].attrs["units"] == "tenths of degrees":
            ds[var].attrs["units"] = "degree"
            attrsbak = ds[var].attrs
            ds[var] = ds[var] / 10
            for k in attrsbak:
                ds[var].attrs[k] = attrsbak[k]

    return ds


def beam2inst(vel, reverse=False, force=False, coord_sys="beam", rotmat=False):
    """
    Rotate velocities from beam to instrument coordinates.

    :param dict adcpo: containing the beam velocity data.
    :param bool reverse: If True, this function performs the inverse rotation (inst->beam).
    :param bool force: When true do not check which coordinate system the data is in prior to performing this rotation.
    """
    if not force:
        if not reverse and coord_sys != "beam":
            raise ValueError("The input must be in beam coordinates.")
        if reverse and coord_sys != "inst":
            raise ValueError("The input must be in inst coordinates.")
    # if 'rotmat' in adcpo['config'].keys():
    #     rotmat = adcpo['config']['rotmat']
    # else:

    cs = "inst"
    if reverse:
        # Can't use transpose because rotation is not between
        # orthogonal coordinate systems
        rotmat = np.linalg.inv(rotmat)
        cs = "beam"
    # raw = adcpo['vel'].transpose()
    raw = np.asmatrix(vel)
    # here I end up with an extra dimension of 4
    # vels = np.einsum('ij,jkl->ikl', rotmat, raw)
    # vels = np.einsum('ij,jk->ik', rotmat, raw)
    vels = np.array(np.asmatrix(rotmat) * raw.transpose())
    # vels = np.einsum('ij,jkl->ikl', rotmat, adcpo['vel'])
    # ValueError: operands could not be broadcast together with remapped
    # shapes [original->remapped]: (4,4)->(4,newaxis,newaxis,4) (16,4,1)->(4,1,16)
    # adcpo['vel'] = np.einsum('ij,jkl->ikl', rotmat, adcpo['vel'])
    return vels.transpose()
    # adcpo['props']['coord_sys'] = cs


def calc_beam_rotmatrix(theta=20, convex=True, degrees=True):
    """
    Calculate the rotation matrix from beam coordinates to
    instrument head coordinates.
    per dolfyn rotate.py code here: https://github.com/lkilcher/dolfyn

    :param float theta: is the angle of the heads (usually 20 or 30 degrees)
    :param int convex: is a flag for convex or concave head configuration.
    :param bool degrees: is a flag which specifies whether theta is in degrees or radians (default: degrees=True)

    """
    deg2rad = np.pi / 180.0
    if degrees:
        theta = theta * deg2rad
    if convex == 0 or convex == -1:
        c = -1
    else:
        c = 1
    a = 1 / (2.0 * np.sin(theta))
    b = 1 / (4.0 * np.cos(theta))
    d = a / (2.0 ** 0.5)
    return np.array(
        [[c * a, -c * a, 0, 0], [0, 0, -c * a, c * a], [b, b, b, b], [d, d, -d, -d]]
    )


def inst2earth(
    vel,
    heading,
    pitch,
    roll,
    reverse=False,
    fixed_orientation=False,
    force=False,
    orientation="up",
):
    """
    Rotate velocities from the instrument to earth coordinates.

    :param dict adcpo: containing the data in instrument coordinates
    :param bool reverse: If True, this function performs the inverse rotation (earth->inst).
    :param bool fixed_orientation: When true, take the average orientation and apply it over the whole record.
    :param bool force: When true do not check which coordinate system the data is in prior to performing this rotation.

    Notes
    -----
    The rotation matrix is taken from the Teledyne RDI ADCP Coordinate Transformation manual January 2008

    When performing the forward rotation, this function sets the 'inst2earth:fixed' flag to the value of
    `fixed_orientation. When performing the reverse rotation, that value is 'popped' from the props dict and the input
    value to this function`fixed_orientation` has no effect. If `'inst2earth:fixed'` is not in the props dict then
    the input value *is* used.
    """
    deg2rad = np.pi / 180.0
    # if not force:
    #     if not reverse and adcpo['props']['coord_sys'] != 'inst':
    #         raise ValueError('The input must be in inst coordinates.')
    #     if reverse and adcpo['props']['coord_sys'] != 'earth':
    #         raise ValueError('The input must be in earth coordinates.')
    # if not reverse and 'declination' in adcpo['props'].keys() and not adcpo['props']['declination_in_heading']:
    #     # Only do this if making the forward rotation.
    #     adcpo['heading_deg'] += adcpo['props']['declination']
    #     adcpo['props']['declination_in_heading'] = True

    r = roll * deg2rad
    p = np.arctan(np.tan(pitch * deg2rad) * np.cos(r))
    h = heading * deg2rad
    # if adcpo['config']['orientation'].lower() == 'up':
    if orientation == "up":
        r += np.pi
    ch = np.cos(h)
    sh = np.sin(h)
    cr = np.cos(r)
    sr = np.sin(r)
    cp = np.cos(p)
    sp = np.sin(p)
    rotmat = np.empty((3, 3, 1))
    rotmat[0, 0, :] = ch * cr + sh * sp * sr
    rotmat[0, 1, :] = sh * cp
    rotmat[0, 2, :] = ch * sr - sh * sp * cr
    rotmat[1, 0, :] = -sh * cr + ch * sp * sr
    rotmat[1, 1, :] = ch * cp
    rotmat[1, 2, :] = -sh * sr - ch * sp * cr
    rotmat[2, 0, :] = -cp * sr
    rotmat[2, 1, :] = sp
    rotmat[2, 2, :] = cp * cr
    # Only operate on the first 3-components, b/c the 4th is err_vel
    # ess = 'ijk,jlk->ilk'
    cs = "earth"
    # if reverse:
    #     cs = 'inst'
    #     fixed_orientation = adcpo['props'].pop('inst2earth:fixed', fixed_orientation)
    #     # ess = ess.replace('ij', 'ji')
    # else:
    #     adcpo['props']['inst2earth:fixed'] = fixed_orientation
    # if fixed_orientation:
    #     # ess = ess.replace('k,', ',')
    #     rotmat = rotmat.mean(-1)
    # todo is the einsum method better?  If so, uncomment the ess statements above
    # vels = np.einsum(ess, rotmat, adcpo['vel'][:,:3])
    vels = np.asmatrix(rotmat) * np.asmatrix(vel[:, :3].transpose())
    vel[:, :3] = vels.transpose()
    # adcpo['props']['coord_sys'] = cs

    return vel
