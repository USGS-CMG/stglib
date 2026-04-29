#!/usr/bin/env python

import matplotlib.pyplot as plt
import pandas as pd

import stglib

args = stglib.cmd.rbrturnaround_parser().parse_args()

df = pd.read_csv(f"{args.basefile}_data.txt", engine="pyarrow")

df = df.rename(columns={"Time": "time"}).set_index("time")

ds = df.to_xarray()

ds.to_netcdf(f"{args.basefile}_turnaround.nc")


for v in ds.data_vars:
    ds[v].plot()
    plt.savefig(f"{args.basefile}_{v}_turnaround.png", dpi=300)

print("Finished creating turnaround plots")
