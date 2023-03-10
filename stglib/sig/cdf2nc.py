import math

import numpy as np
import xarray as xr
from tqdm import tqdm

from ..aqd import aqdutils
from ..core import qaqc, utils


def cdf_to_nc(cdf_filename, atmpres=False):
    """
    Load a "raw" .cdf file and generate a processed .nc file
    """

    # TODO: Add atmospheric pressure offset
    ds = xr.open_dataset(cdf_filename)

    # Clip data to in/out water times or via good_ens
    ds = utils.clip_ds(ds)

    ds = utils.create_nominal_instrument_depth(ds)
    #
    # ds, T, T_orig = set_orientation(ds, ds["TransMatrix"].values)
    #
    e, n, u1, u2 = coord_transform_4beam(
        ds["VelBeam1"].values,
        ds["VelBeam2"].values,
        ds["VelBeam3"].values,
        ds["VelBeam4"].values,
        ds["Heading"].values,
        ds["Pitch"].values,
        ds["Roll"].values,
        ds["Burst_Beam2xyz"].values,
        cs=ds.attrs["SIGBurst_CoordSystem"],
    )

    ds["U"] = xr.DataArray(e, dims=("time", "bindist"))
    ds["V"] = xr.DataArray(n, dims=("time", "bindist"))
    ds["W1"] = xr.DataArray(u1, dims=("time", "bindist"))
    ds["W2"] = xr.DataArray(u2, dims=("time", "bindist"))

    ds = aqdutils.magvar_correct(ds)

    # Rename DataArrays for EPIC compliance
    ds = aqdutils.ds_rename(ds)

    # Add EPIC and CMG attributes
    ds = aqdutils.ds_add_attrs(ds, inst_type="SIG")

    # Add start_time and stop_time attrs
    ds = utils.add_start_stop_time(ds)

    # Add history showing file used
    ds = utils.add_history(ds)

    ds = utils.add_standard_names(ds)

    if "prefix" in ds.attrs:
        nc_filename = ds.attrs["prefix"] + ds.attrs["filename"] + "-a.nc"
    else:
        nc_filename = ds.attrs["filename"] + "-a.nc"

    ds.to_netcdf(nc_filename)
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])

    print("Done writing netCDF file", nc_filename)

    return ds


def coord_transform_4beam(
    vel1, vel2, vel3, vel4, heading, pitch, roll, T, cs="BEAM", out="ENU"
):
    # https://github.com/NortekSupport/coordinatetransforms
    hh = np.pi * (heading - 90) / 180
    pp = np.pi * pitch / 180
    rr = np.pi * roll / 180

    row, col = vel1.shape
    # print(f"{row=}")
    # print(f"{col=}")

    Tmat = np.tile(T, (1, 1, row))

    pmat = aqdutils.make_tilt_np(pp, rr)
    hmat = aqdutils.make_heading_np(hh)

    # print(f"{pmat.shape=}")
    # print(f"{hmat.shape=}")

    r1mat = np.zeros((4, 4, row))

    # print(f"{r1mat.shape=}")

    for i in range(row):
        r1mat[0:3, 0:3, i] = hmat[:, :, i] * pmat[:, :, i]
        r1mat[3, 0:4, i] = r1mat[2, 0:4, i]
        r1mat[0:4, 3, i] = r1mat[0:4, 2, i]

    r1mat[2, 3, :] = 0
    r1mat[3, 2, :] = 0

    rmat = np.zeros_like(r1mat)
    for i in range(row):
        rmat[:, :, i] = r1mat[:, :, i] * Tmat[:, :, i]

    if cs == "BEAM":
        ENU = np.zeros((row, col, 4))

        # print(f"{vel1.shape=}")
        # print('about to print')
        # print(np.array([vel1, vel2, vel3, vel4]).shape)
        # print(f"{row=}")
        # print(f"{col=}")
        # print(f"{rmat.shape=}")
        # print(f"{np.array([vel1[:,0], vel2[:,0], vel3[:,0], vel4[:,0], ]).shape=}")
        # ENU2 = rmat @ np.array([vel1, vel2, vel3, vel4])
        # print(f"{ENU2.shape=}")

        # ENU3 = np.zeros((row, col, 4))
        # for j in range(col):
        #
        #     print(f"{(rmat @ np.array([vel1[:,j], vel2[:,j], vel3[:,j], vel4[:,j], ])).shape=}")
        # print(f"{ENU3.shape=}")

        for i in tqdm(range(row)):
            for j in range(col):
                ENU[i, j, :] = rmat[:, :, i] @ np.array(
                    [
                        vel1[i, j],
                        vel2[i, j],
                        vel3[i, j],
                        vel4[i, j],
                    ]
                )

        # print(f"{ENU.shape=}")
        # print(ENU2.shape)

        # np.testing.assert_array_equal(ENU, ENU2)

        E = ENU[:, :, 0]
        N = ENU[:, :, 1]
        U1 = ENU[:, :, 2]
        U2 = ENU[:, :, 3]

        return E, N, U1, U2
    else:
        raise NotImplementedError(
            f"stglib does not currently support input of {cs} and output of {out} coordinates"
        )
