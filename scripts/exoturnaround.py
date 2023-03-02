#!/usr/bin/env python

import matplotlib.pyplot as plt

import stglib

args = stglib.cmd.exoturnaround_parser().parse_args()

try:
    ds = stglib.exo.read_exo(args.basefile + ".csv", skiprows=args.skiprows)
except UnicodeDecodeError:
    # try reading as Mac OS Western for old versions of Mac Excel
    ds = stglib.exo.read_exo(
        args.basefile + ".csv", skiprows=args.skiprows, encoding="mac-roman"
    )

plt.figure(figsize=(8.5, 11))
plt.subplot(4, 1, 1)
if "Press_dbar" in ds.variables:
    ds["Press_dbar"].plot()
plt.title(args.basefile)

plt.subplot(4, 1, 2)
if "Sal_psu" in ds.variables:
    ds["Sal_psu"].plot()

plt.subplot(4, 1, 3)
if "Temp_°C" in ds.variables:
    ds["Temp_°C"].plot()

plt.subplot(4, 1, 4)
if "Turbidity_NTU" in ds.variables:
    ds["Turbidity_NTU"].plot()
    plt.ylim(0, ds["Turbidity_NTU"].quantile(0.95))

plt.savefig(args.basefile + "_press_sal_temp_turb.png")

if (
    "Chlorophyll_RFU" in ds.variables
    or "BGA-PE_RFU" in ds.variables
    or "ODO_%_sat" in ds.variables
    or "pH" in ds.variables
):
    plt.figure(figsize=(8.5, 11))
    plt.subplot(4, 1, 1)
    if "Chlorophyll_RFU" in ds.variables:
        ds["Chlorophyll_RFU"].plot()
    plt.title(args.basefile)

    plt.subplot(4, 1, 2)
    if "BGA-PE_RFU" in ds.variables:
        ds["BGA-PE_RFU"].plot()

    plt.subplot(4, 1, 3)
    if "ODO_%_sat" in ds.variables:
        ds["ODO_%_sat"].plot()

    plt.subplot(4, 1, 4)
    if "pH" in ds.variables:
        ds["pH"].plot()

    plt.savefig(args.basefile + "_chl_bga_odo_ph.png")

if "Battery_V" in ds.variables or "Wiper_Position_volt" in ds.variables:
    plt.figure(figsize=(8.5, 11))
    plt.subplot(4, 1, 1)
    if "Battery_V" in ds.variables:
        ds["Battery_V"].plot()
    plt.title(args.basefile)

    plt.subplot(4, 1, 2)
    if "Wiper_Position_volt" in ds.variables:
        ds["Wiper_Position_volt"].plot()

    plt.savefig(args.basefile + "_batt_wiper.png")
