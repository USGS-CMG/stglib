import warnings

import pandas as pd
import xarray as xr

from .core import qaqc, utils


def read_exo(filnam, skiprows=8, encoding="utf-8"):
    """Read data from a YSI EXO multiparameter sonde .csv file into an xarray
    Dataset.

    Parameters
    ----------
    filnam : string
        The filename
    skiprows : int, optional
        How many header rows to skip. Default 8
    encoding : string, optional
        File encoding. Default 'utf-8'

    Returns
    -------
    xarray.Dataset
        An xarray Dataset of the EXO data
    """

    try:
        exo = pd.read_csv(filnam, skiprows=skiprows, encoding=encoding, index_col=False)

    # leaving these exceptions, I don't have a mac-roman encoding file to test on though...
    # removed hard coding of parse dates in the exception
    except UnicodeDecodeError:
        exo = pd.read_csv(
            filnam,
            skiprows=skiprows,
            encoding="mac-roman",
        )
    except NotImplementedError as e:
        e.add_note(
            " *** Could not decode file. Try saving the csv file using UTF-8 encoding and retrying"
        )
        raise
    except ValueError as e:
        e.add_note(
            " *** Could not decode header. Have you specified skiprows correctly?"
        )
        raise

    # different software versions may output in all caps or just some caps; convert to all lowercase to make uniform between software versions
    exo.columns = exo.columns.str.lower()

    # find time and date variables to use for converting to datetime; software versions vary in date/time column positions and letter casings
    for var in exo:
        if "time" in var and "fract" not in var:
            time_var = var
        elif "date" in var:
            date_var = var

    # not sure of any other way KOR software outputs date; need to update if this changes in future versions
    if date_var == "date (mm/dd/yyyy)":
        date_format = "%m/%d/%Y"
    else:
        date_format = None

    # if "AM' or 'PM' in time, need to specify with '%I' (0:00-11:59) and '%p' (AM and PM)
    if time_var == "time (hh:mm:ss)" and exo[time_var][0][-1] == "M":
        time_format = "%I:%M:%S %p"
    # if "AM' or 'PM' in time, need to specify with '%I' (0:00-11:59) and '%p' (AM and PM); same with fractions of seconds '%f'
    elif time_var == "time (hh:mm:ss.ttt)" and exo[time_var][0][-1] == "M":
        time_format = "%I:%M:%S.%f %p"
    # if time is 00:00-23:59, specify with '%H'
    elif time_var == "time (hh:mm:ss.ttt)" and exo[time_var][0][-1] != "M":
        time_format = "%H:%M:%S.%f"
    # if time is 00:00-23:59, specify with '%H'
    elif time_var == "time (hh:mm:ss)" and exo[time_var][0][-1] != "M":
        time_format = "%H:%M:%S"
    else:
        time_format = None

    # combine formats so can convert to datetime in one line
    date_time_format = date_format + " " + time_format

    # combine date and time columns and convert to datetime
    exo["time"] = pd.to_datetime(
        exo[date_var] + " " + exo[time_var], format=date_time_format
    )
    exo.set_index("time", inplace=True)

    # specify axis = 1 so columns are droppped; don't need time and date separated
    exo = exo.drop([time_var, date_var], axis=1)

    exo.rename(columns=lambda x: x.replace(" ", "_"), inplace=True)
    exo.rename(columns=lambda x: x.replace("/", "_per_"), inplace=True)

    # New software has time decreasing, need to make sure it's increasing
    exo = xr.Dataset(exo).sortby("time")

    # get sonde and sensor serial numbers from headers
    hdr = read_exo_header(filnam, encoding=encoding, skiprows=skiprows)

    for var in exo:
        if var in hdr:
            exo[var].attrs["sensor_serial_number"] = hdr[var]

    for var in exo:
        if "battery_v" in var:
            sonde_serial_number = exo[var].attrs["sensor_serial_number"]

    if not sonde_serial_number:
        warnings.warn(
            "Can not determine instrument serial number, verify that battery voltage is included in the export file"
        )
        sonde_serial_number = "Not found or specified"

    exo.attrs["serial_number"] = sonde_serial_number
    exo.attrs["instrument_type"] = "YSI EXO2 Multiparameter Sonde"

    if "press_psi_a" in exo:
        pvar_psi = "press_psi_a"
    elif "pressure_psi_a" in exo:
        pvar_psi = "pressure_psi_a"
    else:
        pvar_psi = None
        warnings.warn(
            "*** Could not find pressure (Press_psi_a, Pressure_psi_a) in source data file. Have you exported pressure if this instrument was equipped with a pressure sensor?"
        )

    if pvar_psi:
        # Convert from PSI to dbar
        exo["press_dbar"] = exo[pvar_psi] * 0.689476
        exo["press_dbar"].attrs["sensor_serial_number"] = exo[pvar_psi].attrs[
            "sensor_serial_number"
        ]

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
    Load a raw .cdf file and generate a processed .nc file
    """

    # Load raw .cdf data
    ds = xr.load_dataset(cdf_filename)

    # remove units in case we change and we can use larger time steps
    ds.time.encoding.pop("units")

    # Clip data to in/out water times or via good_ens
    ds = utils.clip_ds(ds)

    ds = ds_rename_vars(ds)

    for k in [
        "press_psi_a",
        "pressure_psi_a",
        "site_name",
        "fault_code",
        "time_(fract._sec)",
        "tds_mg_per_l",
        "tss_mg_per_l",
        "wiper_position_volt",
        "cable_pwr_v",
        # https://www.ysi.com/file%20library/documents/manuals/exo-user-manual-web.pdf
        # nLF_Cond_µS_per_cm: "This convention is typically used in German markets." pp. 85
        "nlf_cond_µs_per_cm",
        "nlf_cond_ms_per_cm",
        "vertical_position_m",
        "ph_mv",
        "file_name",
        "user_id",
        "wiper_position_volt",
        "odo_%_cb",
    ]:
        if k in ds:
            ds = ds.drop(k)

    ds = qaqc.drop_vars(ds)

    if atmpres:
        ds = utils.atmos_correct(ds, atmpres)

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
    if "cond_ms_per_cm" in ds:
        ds["cond_ms_per_cm"].values = (
            ds["cond_ms_per_cm"].values / 10
        )  # convert from mS/cm to S/m
    if "cond_µs_per_cm" in ds:
        ds["cond_µs_per_cm"].values = (
            ds["cond_µs_per_cm"].values / 10000
        )  # convert from µS/cm to S/m
    if "spcond_ms_per_cm" in ds:
        ds["spcond_ms_per_cm"].values = (
            ds["spcond_ms_per_cm"].values / 10
        )  # convert from mS/cm to S/m
    if "spcond_µs_per_cm" in ds:
        ds["spcond_µs_per_cm"].values = (
            ds["spcond_µs_per_cm"].values / 10000
        )  # convert from µS/cm to S/m

    # set up dict of instrument -> EPIC variable names
    varnames = {
        "press_dbar": "P_1",
        "battery_v": "Bat_106",
        "fdom_rfu": "fDOMRFU",
        "fdom_qsu": "fDOMQSU",
        # capitalization based on Chincoteague names
        "chlorophyll_rfu": "CHLrfu",
        "chlorophyll_µg_per_l": "Fch_906",
        "chlorophyll_ug_per_l": "Fch_906",  # added variable name
        "bga-pe_rfu": "TALPErfu",  # BGA is old variable name
        "bga_pe_rfu": "TALPErfu",  # BGA is old variable name
        "bga-pe_µg_per_l": "TALPE",  # BGA is old variable name
        "bga_pe_ug_per_l": "TALPE",  # BGA is old variable name
        "tal_pe_rfu": "TALPErfu",  # added variable name
        "tal_pe_ug_per_l": "TALPE",  # added variable name
        "tal_pe_µg_per_l": "TALPE",  # added variable name; with micro special character
        "temp_°c": "T_28",
        "temp_∞c": "T_28",
        "cond_ms_per_cm": "C_51",
        "cond_µs_per_cm": "C_51",
        "spcond_ms_per_cm": "SpC_48",
        "spcond_µs_per_cm": "SpC_48",
        "sal_psu": "S_41",
        "odo_%_sat": "OST_62",
        "odo_mg_per_l": "DO",
        "turbidity_ntu": "Turb",
        "turbidity_fnu": "Turb_FNU",
        "ph": "pH_159",
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

    ds["Bat_106"].attrs.update(
        {"units": "V", "long_name": "Battery voltage", "epic_code": 106}
    )

    if "fDOMRFU" in ds:
        ds["fDOMRFU"].attrs.update(
            {
                "units": "percent",
                "long_name": "Fluorescent dissolved organic matter, RFU",
                "comment": "Relative fluorescence units (RFU)",
            }
        )

    if "fDOMQSU" in ds:
        ds["fDOMQSU"].attrs.update(
            {
                "units": "1e-9",
                "long_name": "Fluorescent dissolved organic matter, QSU",
                "comment": "Quinine sulfate units (QSU)",
            }
        )

    if "CHLrfu" in ds:
        ds["CHLrfu"].attrs.update(
            {
                "units": "percent",
                "long_name": "Chlorophyll A, RFU",
                "comment": "Relative fluorescence units (RFU)",
            }
        )

    if "Fch_906" in ds:
        ds["Fch_906"].attrs.update(
            {
                "units": "ug/L",
                "long_name": "Chlorophyll A",
                "epic_code": 906,
                "standard_name": "mass_concentration_of_chlorophyll_in_sea_water",
                "comment": "from calibration of sensor with rhodamine W/T in lab",
            }
        )

    if "TALPErfu" in ds:
        ds["TALPErfu"].attrs.update(
            {
                "units": "percent",
                "long_name": "Total algae phycoerythrin, RFU",
                "comment": "Relative fluorescence units (RFU); formerly called BGAPErfu (Blue green algae phycoerythrin, RFU)",
            }
        )

    if "TALPE" in ds:
        ds["TALPE"].attrs.update(
            {
                "units": "ug/L",
                "long_name": "Total algae phycoerythrin",
                "comment": "Formerly called BGAPE (Blue green algae phycoerythrin)",
            }
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
            "units": "S/m",
            "long_name": "Specific Conductivity",
            "comment": "Temperature compensated to 25 °C",
            "epic_code": 48,
            "standard_name": "sea_water_electrical_conductivity_at_reference_temperature",
        }
    )

    ds["S_41"].attrs.update(
        {
            "units": "1",
            "long_name": "Salinity, PSU",
            "comment": "Practical salinity units (PSU)",
            "epic_code": 41,
            "standard_name": "sea_water_practical_salinity",
        }
    )

    if "OST_62" in ds:
        ds["OST_62"].attrs.update(
            {
                "units": "percent",
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
                "comment": "Nephelometric turbidity units (NTU)",
                "standard_name": "sea_water_turbidity",
            }
        )

    if "Turb_FNU" in ds:
        ds["Turb_FNU"].attrs.update(
            {
                "units": "1",
                "long_name": "Turbidity, FNU",
                "comment": "Formazin nephelometric units (FNU)",
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
                "long_name": "Uncorrected pressure",
                "epic_code": 1,
                "standard_name": "sea_water_pressure",
            }
        )

    if "P_1ac" in ds:
        ds["P_1ac"].attrs.update(
            {
                "units": "dbar",
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
                "height_depth_units": "m",
            }
        )

    return ds


def read_exo_header(filnam, skiprows=8, encoding="utf-8"):
    """
    Use header lines initally skipped
    """

    header = {}

    # Old version of KOR export file (headers ~25 lines); updated so it works with latest code version
    with open(filnam, encoding=encoding) as f:
        # make lower case to account for upper/lower case variations in files
        data = f.read().lower()
        if "devices list:" in data:

            data = data.splitlines()

            variables = data[skiprows].split(",")

            # replace spaces and "/" as done with exo pd.DataFrame
            variables = [
                strings.replace(" ", "_").replace("/", "_per_") for strings in variables
            ]

            # copying this and later the variable names will be replaced with serial numbers
            sen_nums = variables.copy()

            # start at -1 one so first row is 0 (python indexing)
            row_count = -1
            for row in data[0:skiprows]:
                row_count += 1
                if "devices list:" in row:
                    # add two more rows; sensors and serial numbers are beginning to be listed three rows after "devices list" row
                    start = row_count + 3

                elif "1,2,3,4,5" in row:
                    # remove one row: sensors and serial numbers stop being listed two rows before this, but indexing below is not inclusive so only remove one
                    end = row_count - 1

            # start is inclusive here, end it not inclusive
            for row in data[start:end]:
                row = row.strip().split(",")
                # if ";" in row[3]:
                if len(row[3]) > 1:
                    cols = row[3].split(";")
                    for k in cols:
                        # subtract one for python indexing
                        sen_nums[int(k) - 1] = row[1]

                elif len(row[3]) == 1:
                    # subtract one for python indexing
                    sen_nums[int(row[3]) - 1] = row[1]

                elif len(row[3]) == 0:
                    # wiper will not have correspond variable (we don't release wiper position data)
                    pass

            # serial numbers are technically upper case, so change to upper
            sen_nums = [strings.upper() for strings in sen_nums]

            header = dict(zip(variables, sen_nums))

        else:
            # Newer versions of KOR export file; headers can have some variations between KOR software versions, this accounts for all the current variations
            data = data.splitlines()
            for row in data:
                if "sensor serial" in row:
                    sen_nums = row.split(",")
                    # serial numbers are technically upper case, so change to upper
                    sen_nums = [strings.upper() for strings in sen_nums]
                if "date" in row:
                    variables = row.split(",")
                    # replace spaces and "/" as done with exo pd.DataFrame
                    variables = [
                        strings.replace(" ", "_").replace("/", "_per_")
                        for strings in variables
                    ]

            header = dict(zip(variables, sen_nums))

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

        ds = qaqc.trim_min_diff_pct(ds, var)

        ds = qaqc.trim_max_diff(ds, var)

        ds = qaqc.trim_max_diff_pct(ds, var)

        ds = qaqc.trim_med_diff(ds, var)

        ds = qaqc.trim_med_diff_pct(ds, var)

        ds = qaqc.trim_bad_ens(ds, var)

    for var in varlist:
        ds = qaqc.trim_by_any(
            ds, var
        )  # re-run and trim by other variables as necessary

    return ds
