import numpy as np
import xarray as xr
from numpy.linalg import inv

from stglib.aqd import aqdutils


def coord_transform(ds, out="enu"):
    """Perform coordinate transformation"""

    cs = ds.attrs["VECCoordinateSystem"]

    if "transform_to" in ds.attrs:
        out = ds.attrs["transform_to"]

    if cs.lower() == out.lower():
        print(
            f"Data already in {cs} coordinates and conversion to {out} requested. Doing nothing."
        )

    elif cs.lower() == "beam" and out.lower() == "xyz":

        # transform to xyz (inst)
        x, y, z = vec_beam_to_xyz(ds.VEL1, ds.VEL2, ds.VEL3, ds["TransMatrix"])

        ds["X"] = xr.DataArray(x, dims=("time"))
        ds["Y"] = xr.DataArray(y, dims=("time"))
        ds["Z"] = xr.DataArray(z, dims=("time"))

    elif cs.lower() == "xyz" and out.lower() == "enu":

        # need to create orienation matrix using heading, pitch, and roll
        orientmat = create_orientmat(ds)

        u, v, w = vec_xyz_to_enu(ds.VEL1, ds.VEL2, ds.VEL3, orientmat)

        ds["U"] = xr.DataArray(u, dims=("time"))
        ds["V"] = xr.DataArray(v, dims=("time"))
        ds["W"] = xr.DataArray(w, dims=("time"))

    elif cs.lower() == "beam" and out.lower() == "enu":

        # first transform to xyz (inst)
        x, y, z = vec_beam_to_xyz(ds.VEL1, ds.VEL2, ds.VEL3, ds["TransMatrix"])

        ds["X"] = xr.DataArray(x, dims=("time"))
        ds["Y"] = xr.DataArray(y, dims=("time"))
        ds["Z"] = xr.DataArray(z, dims=("time"))

        # need to create orienation matrix using heading, pitch, and roll
        orientmat = create_orientmat(ds)

        # now transform to enu
        u, v, w = vec_xyz_to_enu(ds.X, ds.Y, ds.Z, orientmat)

        ds["U"] = xr.DataArray(u, dims=("time"))
        ds["V"] = xr.DataArray(v, dims=("time"))
        ds["W"] = xr.DataArray(w, dims=("time"))

        # do not need intermediate XYZ vels because we keep Beam and now have ENU vels
        ds = ds.drop(["X", "Y", "Z"])

    # reverse transformation enu to xyz
    elif cs.lower() == "enu" and out.lower() == "xyz":

        # need to create orienation matrix using heading, pitch, and roll
        orientmat = create_orientmat(ds)

        # need to transpose orientation matrix for reverse transformation
        orientmat_transpose = orientmat.transpose("inst", "earth", "time")

        # now transform to xyz
        x, y, z = vec_xyz_to_enu(ds.E, ds.N, ds.U, orientmat_transpose)

        ds["X"] = xr.DataArray(x, dims=("time"))
        ds["Y"] = xr.DataArray(y, dims=("time"))
        ds["Z"] = xr.DataArray(z, dims=("time"))

    # reverse transformation enu to beam
    elif cs.lower() == "enu" and out.lower() == "beam":

        # need to create orienation matrix using heading, pitch, and roll
        orientmat = create_orientmat(ds)

        # need to transpose orientation matrix for reverse transformation
        orientmat_transpose = orientmat.transpose("inst", "earth", "time")

        # now transform to xyz
        x, y, z = vec_xyz_to_enu(ds.VEL1, ds.VEL2, ds.VEL3, orientmat_transpose)

        ds["X"] = xr.DataArray(x, dims=("time"))
        ds["Y"] = xr.DataArray(y, dims=("time"))
        ds["Z"] = xr.DataArray(z, dims=("time"))

        T = inv(ds["TransMatrix"])

        # use same def as beam2inst but using inverse, so it's reversed to beam
        v1, v2, v3 = vec_beam_to_xyz(ds.X, ds.Y, ds.Z, T)

        ds["BEAM1VEL"] = xr.DataArray(v1, dims=("time"))
        ds["BEAM2VEL"] = xr.DataArray(v2, dims=("time"))
        ds["BEAM3VEL"] = xr.DataArray(v3, dims=("time"))

    # reserve transformation xyz to beam
    elif cs.lower() == "xyz" and out.lower() == "beam":

        T = inv(ds["TransMatrix"])

        # use same def as beam2inst but using inverse, so it's reversed to beam
        v1, v2, v3 = vec_beam_to_xyz(ds.VEL1, ds.VEL2, ds.VEL3, T)

        ds["BEAM1VEL"] = xr.DataArray(v1, dims=("time"))
        ds["BEAM2VEL"] = xr.DataArray(v2, dims=("time"))
        ds["BEAM3VEL"] = xr.DataArray(v3, dims=("time"))

    return ds


def vec_beam_to_xyz(vel1, vel2, vel3, T):
    vels = np.array([vel1, vel2, vel3])
    xyz = np.einsum("ij,j...->i...", T, vels)

    x = xyz[0]
    y = xyz[1]
    z = xyz[2]

    return x, y, z


def vec_xyz_to_enu(vel1, vel2, vel3, orientmat):
    xyz = np.array([vel1, vel2, vel3])
    uvw = np.einsum("ijk,i...k->j...k", orientmat, xyz)

    u = uvw[0]
    v = uvw[1]
    w = uvw[2]

    return u, v, w


def create_orientmat(ds):

    Heading = ds.Heading.copy()
    Pitch = ds.Pitch.copy()
    Roll = ds.Roll.copy()

    # if bit 0 (starting from right side) in the status code is 1, the instrument was pointing down. Need to change sign of rows 2 and 3 of orientation matrix (orientmat)
    if (ds.orientation).any() == 1:

        Roll[ds.orientation == 1] += 180

    hh = np.pi * (Heading - 90) / 180
    pp = np.pi * Pitch / 180
    rr = np.pi * Roll / 180

    # compute heading and tilt matrix and reorder dimensions so matrix multiplication works properly
    tilt = aqdutils.make_tilt_np(pp, rr)
    tilt = np.moveaxis(tilt, 2, 0)
    head = aqdutils.make_heading_np(hh)
    head = np.moveaxis(head, 2, 0)

    orientmat = head @ tilt
    orientmat = orientmat.transpose(2, 1, 0)

    ds["orientmat"] = xr.DataArray(orientmat, dims=["earth", "inst", "time"])
    ds["orientmat"].attrs[
        "long_name"
    ] = "XYZ (inst) to ENU (earth) Transformation Matrix"
    ds["orientmat"].attrs[
        "note"
    ] = "Generated from instrument heading, pitch, and roll data"

    ds["earth"] = ["E", "N", "U"]
    ds["earth"].attrs["units"] = 1
    ds["earth"].attrs["long_name"] = "Earth Reference Frame"

    return orientmat
