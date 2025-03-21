import operator
import warnings

import numpy as np
import scipy.signal
import xarray as xr

from . import filter, utils


def call_qaqc(ds):
    """Run all QAQC functions. Only those called in the yaml file will be implemented.  If you add a new QAQC function, be sure to add it here too."""

    for var in ds.data_vars:

        ds = filter.apply_butter_filt(ds, var)
        ds = filter.apply_med_filt(ds, var)

        ds = trim_min_diff(ds, var)
        ds = trim_min_diff_pct(ds, var)

        ds = trim_max_diff(ds, var)
        ds = trim_max_diff_pct(ds, var)

        ds = trim_med_diff(ds, var)
        ds = trim_med_diff_pct(ds, var)

        ds = trim_maxabs_diff(ds, var)
        ds = trim_maxabs_diff_2d(ds, var)

        ds = trim_max_blip(ds, var)
        ds = trim_max_blip_pct(ds, var)

        ds = trim_std_ratio(ds, var)
        ds = trim_max_std(ds, var)

        ds = trim_fliers(ds, var)
        ds = trim_warmup(ds, var)

        ds = trim_min(ds, var)
        ds = trim_max(ds, var)

        ds = trim_bad_ens(ds, var)
        ds = trim_bad_ens_indiv(ds, var)

    for var in ds.data_vars:
        ds = trim_mask(ds, var)  # re-run and check for masking
        ds = trim_mask_expr(ds, var)

    for var in ds.data_vars:
        ds = trim_by_any(ds, var)  # re-run and trim by other variables as necessary

    ds = drop_vars(ds)

    return ds


def trim_min(ds, var):
    if var + "_min" in ds.attrs:
        if "sample" in ds:
            cond = (ds[var] >= ds.attrs[var + "_min"]).any(dim="sample")
        else:
            cond = ds[var] >= ds.attrs[var + "_min"]

        affected = cond.size - cond.sum() - ds[var].isnull().sum()
        ds[var] = ds[var].where(cond)

        notetxt = f"Values filled where less than {ds.attrs[var + '_min']} units; {affected.values} values affected. "

        ds = utils.insert_note(ds, var, notetxt)

    return ds


def trim_max(ds, var):
    if var + "_max" in ds.attrs:
        cond = ds[var] <= ds.attrs[var + "_max"]
        affected = cond.size - cond.sum() - ds[var].isnull().sum()
        ds[var] = ds[var].where(cond)

        notetxt = f"Values filled where greater than {ds.attrs[var + '_max']} units; {affected.values} values affected. "

        ds = utils.insert_note(ds, var, notetxt)

    return ds


def trim_min_diff(ds, var):
    if var + "_min_diff" in ds.attrs:
        cond = np.ediff1d(ds[var], to_begin=0) < ds.attrs[var + "_min_diff"]
        affected = cond.sum()
        ds[var][cond] = np.nan

        notetxt = f"Values filled where data decreases by more than {ds.attrs[var + '_min_diff']} units in a single time step; {affected} values affected. "

        ds = utils.insert_note(ds, var, notetxt)

    return ds


def trim_min_diff_pct(ds, var):
    if var + "_min_diff_pct" in ds.attrs:
        cond = (
            100 * np.ediff1d(ds[var], to_begin=0) / np.roll(ds[var], 1)
            < ds.attrs[var + "_min_diff_pct"]
        )
        affected = cond.sum()
        ds[var][cond] = np.nan

        notetxt = f"Values filled where data decreases by more than {ds.attrs[var + '_min_diff_pct']} percent in a single time step; {affected} values affected. "

        ds = utils.insert_note(ds, var, notetxt)

    return ds


def trim_max_diff(ds, var):
    if var + "_max_diff" in ds.attrs:
        cond = np.ediff1d(ds[var], to_begin=0) > ds.attrs[var + "_max_diff"]
        affected = cond.sum()
        ds[var][cond] = np.nan

        notetxt = f"Values filled where data increases by more than {ds.attrs[var + '_max_diff']} units in a single time step; {affected} values affected. "

        ds = utils.insert_note(ds, var, notetxt)

    return ds


def trim_max_diff_pct(ds, var):
    if var + "_max_diff_pct" in ds.attrs:
        cond = (
            100 * np.ediff1d(ds[var], to_begin=0) / np.roll(ds[var], 1)
            > ds.attrs[var + "_max_diff_pct"]
        )
        affected = cond.sum()
        ds[var][cond] = np.nan

        notetxt = f"Values filled where data increases by more than {ds.attrs[var + '_max_diff_pct']} percent in a single time step; {affected} values affected. "

        ds = utils.insert_note(ds, var, notetxt)

    return ds


def trim_med_diff(ds, var):
    if var + "_med_diff" in ds.attrs:
        if "kernel_size" in ds.attrs:
            kernel_size = ds.attrs["kernel_size"]
        else:
            kernel_size = 5
        filtered = scipy.signal.medfilt(ds[var], kernel_size=kernel_size)
        bads = np.abs(ds[var] - filtered) > ds.attrs[var + "_med_diff"]
        affected = bads.sum()
        ds[var][bads] = np.nan

        notetxt = f"Values filled where difference between {kernel_size}-point median filter and original values is greater than {ds.attrs[var + '_med_diff']}; {affected.values} values affected. "

        ds = utils.insert_note(ds, var, notetxt)

    return ds


def trim_med_diff_pct(ds, var):
    if var + "_med_diff_pct" in ds.attrs:
        if "kernel_size" in ds.attrs:
            kernel_size = ds.attrs["kernel_size"]
        else:
            kernel_size = 5
        filtered = scipy.signal.medfilt(ds[var], kernel_size=kernel_size)
        bads = (
            100 * np.abs(ds[var] - filtered) / ds[var] > ds.attrs[var + "_med_diff_pct"]
        )
        affected = bads.sum()
        ds[var][bads] = np.nan

        notetxt = f"Values filled where percent difference between {kernel_size}-point median filter and original values is greater than {ds.attrs[var + '_med_diff_pct']}; {affected.values} values affected. "

        ds = utils.insert_note(ds, var, notetxt)

    return ds


def trim_bad_ens(ds, var):
    if var + "_bad_ens" in ds.attrs:
        inc = np.arange(0, len(ds.attrs[var + "_bad_ens"]), 2)
        for n in inc:
            print(
                "%s: Trimming using bad_ens %s"
                % (var, str(ds.attrs[var + "_bad_ens"][n : n + 2]))
            )
            if isinstance(ds.attrs[var + "_bad_ens"][n], str):
                bads = (ds["time"] >= np.datetime64(ds.attrs[var + "_bad_ens"][n])) & (
                    ds["time"] <= np.datetime64(ds.attrs[var + "_bad_ens"][n + 1])
                )
                ds[var] = ds[var].where(~bads)
            else:
                bads = np.full(ds[var].shape, False)
                bads[
                    np.arange(
                        ds.attrs[var + "_bad_ens"][n], ds.attrs[var + "_bad_ens"][n + 1]
                    )
                ] = True
                ds[var][bads] = np.nan

            notetxt = "Data clipped using bad_ens values of %s. " % (
                str(ds.attrs[var + "_bad_ens"])
            )

            ds = utils.insert_note(ds, var, notetxt)

    return ds


def trim_bad_ens_indiv(ds, var):
    # a list of individual bad points, not a list of ranges
    if var + "_bad_ens_indiv" in ds.attrs:
        trmvars = ds.attrs[var + "_bad_ens_indiv"]

        if not isinstance(trmvars, list):  # make sure it is a list before looping
            trmvars = [trmvars]

        for x in trmvars:
            if isinstance(x, str):
                bads = ds["time"] == np.datetime64(x)
                ds[var] = ds[var].where(~bads)
            else:
                bads = np.full(ds[var].shape, False)
                bads[x] = True
                ds[var] = ds[var].where(~bads)

        notetxt = "Data clipped using bad_ens_indiv values of %s. " % (
            str(ds.attrs[var + "_bad_ens_indiv"])
        )

        ds = utils.insert_note(ds, var, notetxt)

    return ds


def trim_by_any(ds, var):
    attrlist = []
    for a in ds.attrs:
        if "trim_by" in a:
            attrlist.append(a)

    for a in attrlist:
        if a in ds.attrs and ds.attrs[a].lower() == "true" and var in ds:
            # xarray doesn't support writing attributes as booleans
            if f"{a}_exclude" in ds.attrs and var in ds.attrs[f"{a}_exclude"]:
                pass
            else:
                trimvar = a.split("trim_by_")[-1]
                print(f"{var}: Trimming using valid {trimvar} threshold")
                ds[var] = ds[var].where(~ds[trimvar].isnull())

                if var != trimvar:
                    notetxt = f"Values filled using valid {trimvar} threshold. "

                    ds = utils.insert_note(ds, var, notetxt)

    return ds


def trim_max_std(ds, var):
    if var + "_max_std" in ds.attrs:
        if var + "_std" in ds.data_vars:  # check to see that std var exists
            varstd = var + "_std"
            cond = ds[varstd] > ds.attrs[var + "_max_std"]
            affected = cond.sum()
            ds[var][cond] = np.nan

            notetxt = f"Values filled where standard deviation greater than {ds.attrs[var + '_max_std']} units; {affected.values} values affected. "

            ds = utils.insert_note(ds, var, notetxt)

        else:
            print(
                f"{var}_std does not exist was NOT able to trim using maximum standard deviation method"
            )

    return ds


def trim_max_blip(ds, var):
    """trim short-lived maximum "blips", values that increase and then immediately decrease at the next time step"""
    if var + "_max_blip" in ds.attrs:
        cond = (np.ediff1d(ds[var], to_begin=0) > ds.attrs[var + "_max_blip"]) & (
            np.ediff1d(ds[var], to_end=0) < -ds.attrs[var + "_max_blip"]
        )
        affected = cond.sum()
        ds[var][cond] = np.nan

        notetxt = f"Values filled where data blips by more than {ds.attrs[var + '_max_blip']} units in a single time step; {affected} values affected. "

        ds = utils.insert_note(ds, var, notetxt)

    return ds


def trim_max_blip_pct(ds, var):
    """trim short-lived maximum "blips", values that increase and then immediately decrease at the next time step"""
    if var + "_max_blip_pct" in ds.attrs:
        cond = (
            100 * np.ediff1d(ds[var], to_begin=0) / np.roll(ds[var], 1)
            > ds.attrs[var + "_max_blip_pct"]
        ) & (
            100 * np.ediff1d(ds[var], to_end=0) / np.roll(ds[var], -1)
            < -ds.attrs[var + "_max_blip_pct"]
        )
        affected = cond.sum()
        ds[var][cond] = np.nan

        notetxt = f"Values filled where data blips by more than {ds.attrs[var + '_max_blip_pct']} percent in a single time step; {affected} values affected. "

        ds = utils.insert_note(ds, var, notetxt)

    return ds


def trim_fliers(ds, var):
    """trim "fliers", single (or more) presumably bad data points unconnected to other, good data points"""

    if var + "_fliers" in ds.attrs:
        num = ds.attrs[var + "_fliers"]

        a = np.ma.masked_array(ds[var], fill_value=np.nan)
        a[np.isnan(a)] = np.ma.masked
        for b in np.ma.clump_unmasked(a):
            if b.stop - b.start <= num:
                a[b] = np.nan
        a = a.filled()

        cond = np.isfinite(a)
        affected = cond.size - cond.sum() - ds[var].isnull().sum()
        ds[var] = ds[var].where(cond)

        notetxt = f"Fliers of {ds.attrs[var + '_fliers']} or fewer points removed; {affected.values} values affected. "

        ds = utils.insert_note(ds, var, notetxt)

    return ds


def trim_maxabs_diff_2d(ds, var):
    """trim values in a 2D DataArray when the absolute value of the increase is greater than a specified amount"""

    if var + "_maxabs_diff_2d" in ds.attrs:
        try:
            dim1 = ds.attrs[var + "_maxabs_diff_2d"][0]
            val1 = float(ds.attrs[var + "_maxabs_diff_2d"][1])
            dim2 = ds.attrs[var + "_maxabs_diff_2d"][2]
            val2 = float(ds.attrs[var + "_maxabs_diff_2d"][3])

            print(
                f"{var}: Trimming using maximum absolute diff of 2d variable with [dim, val] pairs [{dim1}, {val1}, {dim2}, {val2}]"
            )

            bads1 = np.abs(ds[var].diff(dim=dim1, label="upper")) >= val1
            bads2 = np.abs(ds[var].diff(dim=dim2, label="upper")) >= val2

            # pad bads mask for first element along dim so not set to all NaNs
            padding1 = xr.zeros_like(ds[var].isel({dim1: 0})).astype(bool)
            bads1 = xr.concat([padding1, bads1], dim=dim1)

            padding2 = xr.zeros_like(ds[var].isel({dim2: 0})).astype(bool)
            bads2 = xr.concat([padding2, bads2], dim=dim2)

            ds[var] = ds[var].where(~bads1)
            ds[var] = ds[var].where(~bads2)

            notetxt = f"Values filled where data increases by more than {val1} units (absolute) along {dim1} dim, and {val2} units along {dim2} dim. "
            ds = utils.insert_note(ds, var, notetxt)
        except ValueError:
            raise TypeError(
                f"Values for {var}_maxabs_diff_2d not in required format [dim1(str), val1(float), dim2(str), val2(float)]. No maxabs_diff_2d trimming was done!!"
            )

    return ds


def trim_mask(ds, var):
    """trim values using other variable(s) as mask"""
    if var + "_mask" in ds.attrs:
        trmvars = ds.attrs[var + "_mask"]

        if not isinstance(trmvars, list):  # make sure it is a list before looping
            trmvars = [trmvars]

        for trimvar in trmvars:
            if "beam" in ds[trimvar].dims:
                for bm in ds["beam"].values:
                    cond = ~ds[trimvar].sel(beam=bm).isnull()
                    affected = cond.size - cond.sum() - ds[var].isnull().sum()
                    ds[var] = ds[var].where(cond)

                    notetxt = f"Values filled using {trimvar} beam {bm} mask; {affected.values} values affected. "
                    ds = utils.insert_note(ds, var, notetxt)

            else:
                cond = ~ds[trimvar].isnull()
                affected = cond.size - cond.sum() - ds[var].isnull().sum()
                ds[var] = ds[var].where(cond)

                notetxt = f"Values filled using {trimvar} mask; {affected.values} values affected. "
                ds = utils.insert_note(ds, var, notetxt)

    return ds


def trim_mask_expr(ds, var):
    """
    Trim values based on an expression containing another variable.
    For example, in the config yaml:

    Turb_mask_expr: "P_1ac < 0.1"

    will mask/nan all Turb data where P_1ac is less than 0.1
    """
    operator_map = {
        ">": operator.gt,
        "<": operator.lt,
        ">=": operator.ge,
        "<=": operator.le,
        "==": operator.eq,
        "!=": operator.ne,
    }

    def evaluate_mask_expr(ds, expression):
        parts = expression.split()
        if len(parts) != 3:
            raise ValueError("Invalid mask expression string format")

        left = ds[parts[0]]
        operator_str = parts[1]
        right = float(parts[2])

        if operator_str in operator_map:
            return operator_map[operator_str](left, right)
        else:
            raise ValueError(
                f"Unsupported operator {operator_str}; supported operators are {[x for x in operator_map]}"
            )

    if f"{var}_mask_expr" in ds.attrs:
        mask = ds.attrs[f"{var}_mask_expr"]
        cond = evaluate_mask_expr(ds, mask)
        affected = cond.sum()
        ds[var] = ds[var].where(~cond)
        # print(f"Evaluating '{mask}' with {v} = {ds[v]}: {result}")
        notetxt = (
            f"Values filled using mask of '{mask}'; {affected.values} values affected. "
        )
        ds = utils.insert_note(ds, var, notetxt)

    return ds


def trim_maxabs_diff(ds, var):
    if var + "_maxabs_diff" in ds.attrs:
        if "sample" in ds[var].dims:
            cond = (
                np.abs(ds[var].diff(dim="sample", label="upper"))
                >= ds.attrs[var + "_maxabs_diff"]
            )
            affected = cond.sum()
            val = ds.attrs[var + "_maxabs_diff"]

            # pad bads mask for first element along dim so not set to all NaNs
            padding = xr.zeros_like(ds[var].isel({"sample": 0})).astype(bool)
            cond = xr.concat([padding, cond], dim="sample")
            ds[var] = ds[var].where(~cond)

        else:
            val = ds.attrs[var + "_maxabs_diff"]
            cond = (
                np.abs(np.ediff1d(ds[var], to_begin=0)) > ds.attrs[var + "_maxabs_diff"]
            )
            affected = cond.sum()
            ds[var][cond] = np.nan

        notetxt = f"Values filled where data increases or decreases by more than {val} units in a single time step; {affected} values affected. "

        ds = utils.insert_note(ds, var, notetxt)

    return ds


def trim_std_ratio(ds, var):
    """trim values using ratio of standard deviation value to variable value (requires var_std to exists)"""
    if var + "_std_ratio" in ds.attrs:
        if var + "_std" in ds.data_vars:  # check to see that std var exists
            print(
                "%s: Trimming using standard deviation ratio of %f"
                % (var, ds.attrs[var + "_std_ratio"])
            )

            varstd = var + "_std"
            cond = ds[varstd] / ds[var] > ds.attrs[var + "_std_ratio"]
            affected = cond.sum()
            ds[var][cond] = np.nan

            notetxt = f"Values filled where standard deviation ratio threshold of {ds.attrs[var + '_std_ratio']} was exceeded; {affected.values} values affected. "

            ds = utils.insert_note(ds, var, notetxt)

        else:
            raise ValueError(
                f"User specified {ds.attrs[var + '_std_ratio']=} but {var}_std does not exist. Was not able to trim using standard deviation ratio method"
            )

    return ds


def trim_warmup(ds, var):
    if var + "_warmup_samples" in ds.attrs:
        if "sample" in ds[var].coords:
            ds[var] = ds[var].where(ds["sample"] > ds.attrs[var + "_warmup_samples"])
            notetxt = f"Removed {ds.attrs[var + '_warmup_samples']} samples at the beginning of each burst. "

            ds = utils.insert_note(ds, var, notetxt)
        else:
            raise ValueError(
                f"User specified {ds.attrs[var + '_warmup_samples']=} but {var=} does not have coordinates that include samples."
            )

    return ds


def drop_vars(ds):
    """Remove variables in the final Dataset as specified by the user"""
    dropped = []
    if "drop_vars" in ds.attrs:
        drpvars = ds.attrs["drop_vars"]
        if not isinstance(drpvars, list):  # make sure it is a list before looping
            drpvars = [drpvars]
        for k in drpvars:
            if k in ds:
                ds = ds.drop_vars(k)
                dropped.append(k)
            else:
                warnings.warn(f"{k} not in dataset, cannot drop")
        if dropped:
            print(f"Dropped {dropped} from dataset at user request")

    return ds
