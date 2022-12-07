import numpy as np
import pandas as pd
import xarray as xr

from .core import qaqc, utils


def read_exo(filnam, skiprows=25, encoding="utf-8"):
    """Read data from a YSI EXO multiparameter sonde .csv file into an xarray
    Dataset.

    Parameters
    ----------
    filnam : string
        The filename
    skiprows : int, optional
        How many header rows to skip. Default 25
    encoding : string, optional
        File encoding. Default 'utf-8'

    Returns
    -------
    xarray.Dataset
        An xarray Dataset of the EXO data
    """

    try:
        exo = pd.read_csv(
            filnam,
            skiprows=skiprows,
            infer_datetime_format=True,
            # parse_dates=[['Date (MM/DD/YYYY)',
            #               'Time (HH:MM:SS)']],
            parse_dates=[[0, 1]],
            encoding=encoding,
        )
    except UnicodeDecodeError:
        exo = pd.read_csv(
            filnam,
            skiprows=skiprows,
            infer_datetime_format=True,
            # parse_dates=[['Date (MM/DD/YYYY)',
            #               'Time (HH:MM:SS)']],
            parse_dates=[[0, 1]],
            encoding="mac-roman",
        )
    except NotImplementedError as e:
        print(
            (
                " *** Could not decode file. Try saving the csv file using "
                "UTF-8 encoding and retrying\n"
            ),
            e,
        )
    except ValueError as e:
        print(
            (
                " *** Could not decode header. "
                "Have you specified skiprows correctly?\n"
            ),
            e,
        )
    # exo.rename(columns={'Date (MM/DD/YYYY)_Time (HH:MM:SS)': 'time'},
    #            inplace=True)
    exo.rename(
        columns={exo.columns[0]: "time"}, inplace=True
    )  # rename first column to time.
    # Need to do this because the format of the date/time header can change between versions
    exo.set_index("time", inplace=True)
    exo.rename(columns=lambda x: x.replace(" ", "_"), inplace=True)
    exo.rename(columns=lambda x: x.replace("/", "_per_"), inplace=True)
    pvar = None
    if "Press_psi_a" in exo.columns:
        pvar = "Press_psi_a"
    elif "Pressure_psi_a" in exo.columns:
        pvar = "Pressure_psi_a"
    else:
        print(
            "*** Could not find pressure (Press_psi_a, Pressure_psi_a) in source data file. Have you exported pressure if this instrument was equipped with a pressure sensor?"
        )
    if pvar:
        exo["Press_dbar"] = exo[pvar] * 0.689476

    exo = xr.Dataset(exo)
    hdr = read_exo_header(filnam, encoding=encoding)
    exo.attrs["serial_number"] = hdr["serial_number"]
    exo.attrs["INST_TYPE"] = "YSI EXO2 Multiparameter Sonde"
    exo.attrs["COMPOSITE"] = np.int32(0)

    # Apply sensor serial numbers to each sensor
    for k in exo.variables:
        if "fDOM" in k:
            if "fDOM" in hdr:
                hdrvar = "fDOM"
            elif "fDOM QSU" in hdr:
                hdrvar = "fDOM QSU"
            exo[k].attrs["sensor_serial_number"] = hdr[hdrvar]["sensor_serial_number"]
        elif "Chlorophyll" in k or "BGA" in k or "TAL" in k:
            if "Total Algae BGA-PE" in hdr:
                hdrvar = "Total Algae BGA-PE"
            elif "BGA PE RFU" in hdr:
                hdrvar = "BGA PE RFU"
            elif "TAL PE RFU" in hdr:
                hdrvar = "TAL PE RFU"
            exo[k].attrs["sensor_serial_number"] = hdr[hdrvar]["sensor_serial_number"]
        elif "Temp" in k or "Cond" in k or "Sal" in k:
            if "Unknown CT" in hdr:
                exo[k].attrs["sensor_serial_number"] = hdr["Unknown CT"][
                    "sensor_serial_number"
                ]
            elif "Wiped CT" in hdr:
                exo[k].attrs["sensor_serial_number"] = hdr["Wiped CT"][
                    "sensor_serial_number"
                ]
            else:
                hdrvar = "Sal psu"
            exo[k].attrs["sensor_serial_number"] = hdr[hdrvar]["sensor_serial_number"]
        elif "ODO" in k:
            try:
                exo[k].attrs["sensor_serial_number"] = hdr["Optical DO"][
                    "sensor_serial_number"
                ]
            except KeyError:
                exo[k].attrs["sensor_serial_number"] = hdr["ODO % sat"][
                    "sensor_serial_number"
                ]
        elif k == "Turbidity":
            exo[k].attrs["sensor_serial_number"] = hdr["Turbidity"][
                "sensor_serial_number"
            ]
        elif k == "Turbidity_NTU":
            exo[k].attrs["sensor_serial_number"] = hdr["Turbidity NTU"][
                "sensor_serial_number"
            ]
        elif k == "Turbidity_FNU":
            exo[k].attrs["sensor_serial_number"] = hdr["Turbidity FNU"][
                "sensor_serial_number"
            ]
        elif "pH" in k:
            exo[k].attrs["sensor_serial_number"] = hdr["pH"]["sensor_serial_number"]
        elif "Press" in k or "Depth" in k:
            if "Depth Non-Vented 0-10m" in hdr:
                hdrvar = "Depth Non-Vented 0-10m"
            elif "Depth m" in hdr:
                hdrvar = "Depth m"
            elif "Pressure psi a" in hdr:
                hdrvar = "Pressure psi a"
            else:
                hdrvar = None
            exo[k].attrs["sensor_serial_number"] = hdr[hdrvar]["sensor_serial_number"]

    return exo


def csv_to_cdf(metadata):
    """
    Process EXO .csv file to a raw .cdf file
    """

    basefile = metadata["basefile"]

    try:
        ds = read_exo(basefile + ".csv", skiprows=metadata["skiprows"])
    except UnicodeDecodeError:
        # try reading as Mac OS Western for old versions of Mac Excel
        ds = read_exo(
            basefile + ".csv", skiprows=metadata["skiprows"], encoding="mac-roman"
        )

    metadata.pop("skiprows")

    # write out metadata first, then deal exclusively with xarray attrs
    ds = utils.write_metadata(ds, metadata)

    del metadata

    ds = utils.ensure_cf(ds)

    ds = utils.shift_time(ds, 0)

    # configure file
    cdf_filename = ds.attrs["filename"] + "-raw.cdf"

    ds.to_netcdf(cdf_filename, unlimited_dims=["time"])

    print("Finished writing data to %s" % cdf_filename)

    return ds


def cdf_to_nc(cdf_filename, atmpres=False):
    """
    Load a "raw" .cdf file and generate a processed .nc file
    """

    # Load raw .cdf data
    ds = xr.load_dataset(cdf_filename)
    ds.time.encoding.pop(
        "units"
    )  # remove units in case we change and we can use larger time steps

    # Clip data to in/out water times or via good_ens
    ds = utils.clip_ds(ds)

    ds = ds_rename_vars(ds)

    # ds = ds_add_attrs(ds)

    for k in [
        "Press_psi_a",
        "Pressure_psi_a",
        "Site_Name",
        "Fault_Code",
        "Time_(Fract._Sec)",
        "TDS_mg_per_L",
        "TSS_mg_per_L",
        "Wiper_Position_volt",
        "Cable_Pwr_V",
        # https://www.ysi.com/file%20library/documents/manuals/exo-user-manual-web.pdf
        # nLF_Cond_µS_per_cm: "This convention is typically used in German markets." pp. 85
        "nLF_Cond_µS_per_cm",
        "Vertical_Position_m",
        "pH_mV",
    ]:
        if k in ds:
            ds = ds.drop(k)

    if "drop_vars" in ds.attrs:
        for k in ds.attrs["drop_vars"]:
            if k in ds:
                ds = ds.drop(k)

    if atmpres:
        print("Atmospherically correcting data")

        met = xr.load_dataset(atmpres)
        # need to save attrs before the subtraction, otherwise they are lost
        attrs = ds["P_1"].attrs
        ds["P_1ac"] = ds["P_1"] - met["atmpres"] - met["atmpres"].offset
        print("Correcting using offset of %f" % met["atmpres"].offset)
        ds["P_1ac"].attrs = attrs

    ds = exo_qaqc(ds)

    # assign min/max:
    ds = utils.add_min_max(ds)

    ds = utils.add_start_stop_time(ds)

    ds = utils.add_delta_t(ds)

    # add lat/lon coordinates
    ds = utils.ds_add_lat_lon(ds)

    # ds = utils.create_water_depth(ds)
    ds = utils.create_nominal_instrument_depth(ds)

    ds = utils.create_z(ds)

    ds = ds_add_attrs(ds)

    # add lat/lon coordinates to each variable
    # for var in ds.variables:
    #     if (var not in ds.coords) and ("time" not in var):
    #         # ds = utils.add_lat_lon(ds, var)
    #         # ds = utils.no_p_add_depth(ds, var)
    #         ds = utils.add_z_if_no_pressure(ds, var)
    #         # cast as float32
    #         # ds = utils.set_var_dtype(ds, var)

    # No longer report depth
    if "Depth_m" in ds:
        ds = ds.drop("Depth_m")

    # Write to .nc file
    print("Writing cleaned/trimmed data to .nc file")
    nc_filename = ds.attrs["filename"] + "-a.nc"

    ds.to_netcdf(
        nc_filename, unlimited_dims=["time"], encoding={"time": {"dtype": "i4"}}
    )
    utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])
    print("Done writing netCDF file", nc_filename)


def ds_rename_vars(ds):

    if "Cond_mS_per_cm" in ds:
        ds["Cond_mS_per_cm"].values = (
            ds["Cond_mS_per_cm"].values / 10
        )  # convert from mS/cm to S/m
    if "Cond_µS_per_cm" in ds:
        ds["Cond_µS_per_cm"].values = (
            ds["Cond_µS_per_cm"].values / 10000
        )  # convert from µS/cm to S/m
    if "SpCond_µS_per_cm" in ds:
        ds["SpCond_µS_per_cm"].values = (
            ds["SpCond_µS_per_cm"].values / 10
        )  # convert from µS/cm to mS/cm

    # set up dict of instrument -> EPIC variable names
    varnames = {
        "Press_dbar": "P_1",
        "Battery_V": "Bat_106",
        "fDOM_RFU": "fDOMRFU",
        "fDOM_QSU": "fDOMQSU",
        # capitalization based on Chincoteague names
        "Chlorophyll_RFU": "CHLrfu",
        "Chlorophyll_µg_per_L": "Fch_906",
        "Chlorophyll_ug_per_L": "Fch_906",  # added variable name
        "BGA-PE_RFU": "BGAPErfu",
        "BGA_PE_RFU": "BGAPErfu",  # added variable name
        "BGA-PE_µg_per_L": "BGAPE",
        "BGA_PE_ug_per_L": "BGAPE",  # added variable name
        "TAL_PE_RFU": "TALPErfu",  # added variable name
        "TAL_PE_ug_per_L": "TALPE",  # added variable name
        "Temp_°C": "T_28",
        "Temp_∞C": "T_28",
        "Cond_mS_per_cm": "C_51",
        "Cond_µS_per_cm": "C_51",
        "SpCond_mS_per_cm": "SpC_48",
        "SpCond_µS_per_cm": "SpC_48",
        "Sal_psu": "S_41",
        "ODO_%_sat": "OST_62",
        "ODO_mg_per_L": "DO",
        "Turbidity_NTU": "Turb",
        "Turbidity_FNU": "Turb_FNU",
        "pH": "pH_159",
    }

    # check to make sure they exist before trying to rename
    newvars = {}
    for k in varnames:
        if k in ds:
            newvars[k] = varnames[k]

    return ds.rename(newvars)


def ds_add_attrs(ds):
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

    ds["Bat_106"].attrs.update(
        {"units": "V", "long_name": "Battery voltage", "epic_code": 106}
    )

    if "fDOMRFU" in ds:
        ds["fDOMRFU"].attrs.update(
            {
                "units": "percent",
                "long_name": "Fluorescent dissolved organic matter, RFU",
                "comments": "Relative fluorescence units (RFU)",
            }
        )

    if "fDOMQSU" in ds:
        ds["fDOMQSU"].attrs.update(
            {
                "units": "1e-9",
                "long_name": "Fluorescent dissolved organic matter, QSU",
                "comments": "Quinine sulfate units (QSU)",
            }
        )

    if "CHLrfu" in ds:
        ds["CHLrfu"].attrs.update(
            {
                "units": "percent",
                "long_name": "Chlorophyll A, RFU",
                "comments": "Relative fluorescence units (RFU)",
            }
        )

    if "Fch_906" in ds:
        ds["Fch_906"].attrs.update(
            {"units": "ug/L", "long_name": "Chlorophyll A", "epic_code": 906}
        )

    if "BGAPErfu" in ds:
        ds["BGAPErfu"].attrs.update(
            {
                "units": "percent",
                "long_name": "Blue green algae phycoerythrin, RFU",
                "comments": "Relative fluorescence units (RFU)",
            }
        )

    if "BGAPE" in ds:
        ds["BGAPE"].attrs.update(
            {"units": "ug/L", "long_name": "Blue green algae phycoerythrin"}
        )

    if "TALPErfu" in ds:
        ds["TALPErfu"].attrs.update(
            {
                "units": "percent",
                "long_name": "Total algae phycoerythrin, RFU",
                "comments": "Relative fluorescence units (RFU)",
            }
        )

    if "TALPE" in ds:
        ds["TALPE"].attrs.update(
            {"units": "ug/L", "long_name": "Total algae phycoerythrin"}
        )

    ds["T_28"].attrs.update(
        {
            "units": "degree_C",
            "long_name": "Temperature",
            "epic_code": 28,
            "standard_name": "sea_water_temperature",
        }
    )

    ds["C_51"].attrs.update(
        {
            "units": "S/m",
            "long_name": "Conductivity",
            "epic_code": 51,
            "standard_name": "sea_water_electrical_conductivity",
        }
    )

    ds["SpC_48"].attrs.update(
        {
            "units": "mS/cm",
            "long_name": "Conductivity",
            "comment": "Temperature compensated to 25 °C",
            "epic_code": 48,
            "standard_name": "sea_water_electrical_conductivity",
        }
    )

    ds["S_41"].attrs.update(
        {
            "units": "1e-3",
            "long_name": "Salinity, PSU",
            "comments": "Practical salinity units (PSU)",
            "epic_code": 41,
            "standard_name": "sea_water_salinity",
        }
    )

    if "OST_62" in ds:
        ds["OST_62"].attrs.update(
            {
                "units": "%",
                "long_name": "Oxygen percent saturation",
                "epic_code": 62,
                "standard_name": "fractional_saturation_of_oxygen_in_sea_water",
            }
        )

    if "DO" in ds:
        ds["DO"].attrs.update(
            {
                "units": "mg/L",
                "long_name": "Dissolved oxygen",
                "standard_name": "mass_concentration_of_oxygen_in_sea_water",
            }
        )

    if "Turb" in ds:
        ds["Turb"].attrs.update(
            {
                "units": "1",
                "long_name": "Turbidity, NTU",
                "comments": "Nephelometric turbidity units (NTU)",
                "standard_name": "sea_water_turbidity",
            }
        )

    if "Turb_FNU" in ds:
        ds["Turb_FNU"].attrs.update(
            {
                "units": "1",
                "long_name": "Turbidity, FNU",
                "comments": "Formazin nephelometric units (FNU)",
                "standard_name": "sea_water_turbidity",
            }
        )

    if "pH_159" in ds.variables:
        ds["pH_159"].attrs.update(
            {
                "units": "1",
                "standard_name": "sea_water_ph_reported_on_total_scale",
                "epic_code": 159,
            }
        )

    if "P_1" in ds:
        ds["P_1"].attrs.update(
            {
                "units": "dbar",
                "long_name": "Pressure",
                "epic_code": 1,
                "standard_name": "sea_water_pressure",
            }
        )

    if "P_1ac" in ds:
        ds["P_1ac"].attrs.update(
            {
                "units": "dbar",
                "name": "Pac",
                "long_name": "Corrected pressure",
                "standard_name": "sea_water_pressure_due_to_sea_water",
            }
        )
        if "P_1ac_note" in ds.attrs:
            ds = utils.insert_note(ds, "P_1ac", ds.attrs["P_1ac_note"] + " ")

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

    return ds


def read_exo_header(filnam, encoding="utf-8"):
    header = {}
    try:
        # Old version of KOR export file
        hdr = pd.read_csv(filnam, skiprows=None, encoding=encoding)
        hdr = pd.DataFrame(hdr.iloc[:, 0:4])
        # print(hdr)
        header["serial_number"] = (
            hdr[hdr["KOR Export File"] == "Sonde ID"].values[0][1].split(" ")[1]
        )
        for var in [
            "fDOM",
            "Total Algae BGA-PE",
            "Wiped CT",
            "Unknown CT",
            "Optical DO",
            "Turbidity",
            "pH",
            "Depth Non-Vented 0-10m",
        ]:
            vals = hdr[hdr["KOR Export File"] == var]
            if not vals.empty:
                header[var] = {}
                header[var]["sensor_serial_number"] = vals.values[0][1]
                header[var]["data_columns"] = [
                    int(x) for x in vals.values[0][3].split(";")
                ]
    except (pd.errors.ParserError, KeyError):
        # new version of KOR export file
        hdr = pd.read_csv(filnam, skiprows=4, encoding=encoding)
        hdr = pd.DataFrame(hdr.iloc[:, 3:-1])
        # get instrument SN from filename -- this will fail on files not named according to the Kor default file-naming convention but I'm not sure how else to get it
        header["serial_number"] = filnam.split("/")[-1].split("_")[1]
        row = np.where(hdr.iloc[:, 0] == "SENSOR SERIAL NUMBER:")
        a = np.vstack([hdr.iloc[row[0] + 1, :].values, hdr.iloc[row[0], :].values]).T
        for v in a:
            if v[0] != "Site Name":
                header[v[0]] = {}
                header[v[0]]["sensor_serial_number"] = v[1]

    return header


def exo_qaqc(ds):
    """
    QA/QC
    Trim EXO data based on metadata
    """

    varlist = [
        "S_41",
        "C_51",
        "SpC_48",
        "T_28",
        "Turb",
        "fDOMRFU",
        "fDOMQSU",
        "CHLrfu",
        "Fch_906",
        "BGAPErfu",
        "BGAPE",
        "TALPErfu",
        "TALPE",
        "OST_62",
        "DO",
        "pH_159",
        "P_1ac",
        "P_1",
    ]

    [varlist.append(k) for k in ds.data_vars if k not in varlist]

    for var in varlist:
        ds = qaqc.trim_min(ds, var)

        ds = qaqc.trim_max(ds, var)

        ds = qaqc.trim_min_diff(ds, var)

        ds = qaqc.trim_max_diff(ds, var)

        ds = qaqc.trim_med_diff(ds, var)

        ds = qaqc.trim_med_diff_pct(ds, var)

        ds = qaqc.trim_bad_ens(ds, var)

    for var in varlist:
        ds = qaqc.trim_by_any(
            ds, var
        )  # re-run and trim by other variables as necessary

    return ds
