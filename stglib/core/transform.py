import numpy as np
import xarray as xr

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

        ds["U"] = ds.VEL1
        ds["V"] = ds.VEL2
        ds["W"] = ds.VEL3

    elif cs.lower() == "beam" and out.lower() == "xyz":

        x, y, z = vec_beam_to_xyz(
            ds.VEL1, ds.VEL2, ds.VEL3, ds["TransMatrix"]
        )  # transform to xyz (inst)

        ds["X"] = xr.DataArray(x, dims=("time"))
        ds["Y"] = xr.DataArray(y, dims=("time"))
        ds["Z"] = xr.DataArray(z, dims=("time"))

    elif cs.lower() == "xyz" and out.lower() == "enu":

        orientmat = create_orientmat(
            ds
        )  # need to create orienation matrix using heading, pitch, and roll

        u, v, w = vec_xyz_to_enu(ds.VEL1, ds.VEL2, ds.VEL3, orientmat)

        ds["U"] = xr.DataArray(u, dims=("time"))
        ds["V"] = xr.DataArray(v, dims=("time"))
        ds["W"] = xr.DataArray(w, dims=("time"))

    elif cs.lower() == "beam" and out.lower() == "enu":

        x, y, z = vec_beam_to_xyz(
            ds.VEL1, ds.VEL2, ds.VEL3, ds["TransMatrix"]
        )  # first transform to xyz (inst)

        ds["X"] = xr.DataArray(x, dims=("time"))
        ds["Y"] = xr.DataArray(y, dims=("time"))
        ds["Z"] = xr.DataArray(z, dims=("time"))

        orientmat = create_orientmat(
            ds
        )  # need to create orienation matrix using heading, pitch, and roll

        u, v, w = vec_xyz_to_enu(ds.X, ds.Y, ds.Z, orientmat)  # now transform to enu

        ds["U"] = xr.DataArray(u, dims=("time"))
        ds["V"] = xr.DataArray(v, dims=("time"))
        ds["W"] = xr.DataArray(w, dims=("time"))

    elif (
        cs.lower() == "enu" and out.lower() == "xyz"
    ):  # reverse transformation enu to xyz

        orientmat = create_orientmat(
            ds
        )  # need to create orienation matrix using heading, pitch, and roll

        orientmat_transpose = orientmat.transpose(
            "inst", "earth", "time"
        )  # need to transpose orientation matrix for reverse transformation

        x, y, z = vec_xyz_to_enu(
            ds.E, ds.N, ds.U, orientmat_transpose
        )  # now transform to xyz

        ds["X"] = xr.DataArray(x, dims=("time"))
        ds["Y"] = xr.DataArray(y, dims=("time"))
        ds["Z"] = xr.DataArray(z, dims=("time"))

    elif (
        cs.lower() == "enu" and out.lower() == "beam"
    ):  # reverse transformation enu to beam

        orientmat = create_orientmat(
            ds
        )  # need to create orienation matrix using heading, pitch, and roll

        orientmat_transpose = orientmat.transpose(
            "inst", "earth", "time"
        )  # need to transpose orientation matrix for reverse transformation

        x, y, z = vec_xyz_to_enu(
            ds.VEL1, ds.VEL2, ds.VEL3, orientmat_transpose
        )  # now transform to xyz

        ds["X"] = xr.DataArray(x, dims=("time"))
        ds["Y"] = xr.DataArray(y, dims=("time"))
        ds["Z"] = xr.DataArray(z, dims=("time"))

        T = inv(ds["TransMatrix"])

        v1, v2, v3 = vec_beam_to_xyz(
            ds.X, ds.Y, ds.Z, T
        )  # use same def as beam2inst but using inverse, so it's reversed to beam

        ds["BEAM1VEL"] = xr.DataArray(v1, dims=("time"))
        ds["BEAM2VEL"] = xr.DataArray(v2, dims=("time"))
        ds["BEAM3VEL"] = xr.DataArray(v3, dims=("time"))

    elif (
        cs.lower() == "xyz" and out.lower() == "beam"
    ):  # reserve transformation xyz to beam

        T = inv(ds["TransMatrix"])

        v1, v2, v3 = vec_beam_to_xyz(
            ds.VEL1, ds.VEL2, ds.VEL3, T
        )  # use same def as beam2inst but using inverse, so it's reversed to beam

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

    if (
        ds.orientation
    ).any() == 1:  # if bit 0 (starting from right side) in the status code is 1, the instrument was pointing down. Need to change sign of rows 2 and 3 of orientation matrix (orientmat)
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

    ds["earth"] = ["E", "N", "U"]
    ds["earth"].attrs["units"] = 1
    ds["earth"].attrs["long_name"] = "Earth Reference Frame"

    ds["inst"] = ["X", "Y", "Z"]
    ds["inst"].attrs["units"] = 1
    ds["inst"].attrs["long_name"] = "Inst Reference Frame"

    return orientmat
