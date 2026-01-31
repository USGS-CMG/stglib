from . import utils


def ds_add_attrs(ds):
    """
    Add attributes: units, standard name from CF website, long names
    """

    ds = utils.ds_coord_no_fillvalue(ds)

    ds["time"].attrs.update(
        {"standard_name": "time", "axis": "T", "long_name": "time (UTC)"}
    )

    if "abs" in ds:
        ds["abs"].attrs.update(
            {
                "units": "normalized counts",
                "long_name": "Acoustic backscatter strength",
                "transducer_offset_from_bottom": ds.attrs["initial_instrument_height"],
            }
        )

    if "AGC_1202" in ds:
        ds["AGC_1202"].attrs.update(
            {
                "units": "counts",
                "long_name": "Average Echo Intensity",
            }
        )

    # Instrument type abss
    if "amp" in ds:
        ds["amp"].attrs.update(
            {
                "units": "decibels",
                "long_name": "Acoustic signal amplitude",
                "standard_name": "sound_intensity_level_in_water",
                "transducer_offset_from_bottom": ds.attrs["initial_instrument_height"],
                "note": "abs data converted from counts to decibels using equation: decibels = 20*log10(counts).",
            }
        )

    # Instrument type vec
    if "amp" in ds:
        ds["amp"].attrs.update(
            {
                "units": "Counts",
                "standard_name": "signal_intensity_from_multibeam_acoustic_doppler_velocity_sensor_in_sea_water",
                "long_name": "Acoustic Signal Amplitude",
            }
        )

    # Instrument type vec
    if "amp_avg" in ds:
        ds["amp_avg"].attrs.update(
            {
                "units": "Counts",
                "standard_name": "signal_intensity_from_multibeam_acoustic_doppler_velocity_sensor_in_sea_water",
                "long_name": "Average Acoustic Signal Amplitude",
            }
        )

    if "AMP_723" in ds:
        ds["AMP_723"].attrs.update(
            {
                "units": "percent",
                "long_name": "Acoustic Signal Amplitude Strength",
            }
        )

    if "Bat_106" in ds:
        ds["Bat_106"].attrs.update({"units": "V", "long_name": "Battery voltage"})

    if "bin_depth" in ds and ds.get("instrument_type", "").lower() == "aquascat1000r":
        ds["bin_depth"].attrs.update(
            {
                "units": "m",
                "long_name": "bin depth",
                "bin_size": ds.attrs["ABSAbsBinLengthMM"][0] * 0.001,
                "bin_count": ds.attrs["ABSAbsNumBins"][0],
            }
        )

    if "bindist" in ds and ds.get("instrument_type", "").lower() == "aquascat1000r":
        ds["bindist"].attrs.update(
            {
                "units": "m",
                "long_name": "distance from transducer head",
                "bin_size": ds.attrs["ABSAbsBinLengthMM"][0] * 0.001,
                "bin_count": ds.attrs["ABSAbsNumBins"][0],
            }
        )

    if "BPR_915" in ds:
        ds["BPR_915"].attrs.update(
            {"units": "pascals", "standard_name": "air_pressure"}
        )

    if "brange" in ds:
        ds["brange"].attrs.update(
            {
                "units": "m",
                "long_name": "sensor range to boundary",
                "standard_name": "altimeter_range",
            }
        )

    if "C_51" in ds:
        ds["C_51"].attrs.update(
            {
                "units": "S/m",
                "long_name": "Conductivity",
                "standard_name": "sea_water_electrical_conductivity",
            }
        )

    if "CD_310" in ds:
        ds["CD_310"].attrs.update(
            {
                "units": "degree",
                "long_name": "Current Direction (True)",
                "standard_name": "sea_water_velocity_to_direction",
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

    if "CS_300" in ds:
        ds["CS_300"].attrs.update(
            {
                "units": "m s^-1",
                "long_name": "Current Speed",
                "standard_name": "sea_water_speed",
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

    if "Fch_906" in ds:
        ds["Fch_906"].attrs.update(
            {
                "units": "ug/L",
                "long_name": "Chlorophyll A",
                "standard_name": "mass_concentration_of_chlorophyll_in_sea_water",
                "comment": "from calibration of sensor with rhodamine W/T in lab",
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

    if "fDOMRFU" in ds:
        ds["fDOMRFU"].attrs.update(
            {
                "units": "percent",
                "long_name": "Fluorescent dissolved organic matter, RFU",
                "comment": "Relative fluorescence units (RFU)",
            }
        )

    if "GyroHeading" in ds:
        ds["GyroHeading"].attrs.update(
            {
                "units": "degrees",
                "long_name": "Gyro Heading",
                "standard_name": "platform_orientation",
            }
        )

    if "Hdg_1215" in ds:
        ds["Hdg_1215"].attrs.update(
            {
                "units": "degrees",
                "long_name": "Instrument Heading",
                "note": "Heading is in degrees true. Converted from degrees magnetic using magnetic variation.",
                "standard_name": "platform_orientation",
            }
        )

    if "HeadAngle" in ds:
        ds["HeadAngle"].attrs.update(
            {
                "units": "degrees",
                "long_name": "Transducer head angle",
                "comment": "Angle = 0.3 x (Head Position - 600)",
            }
        )

    if "HeadPosition" in ds:
        ds["HeadPosition"].attrs.update(
            {
                "units": "1",
                "long_name": "Transducer head position",
                "comment": "Integer values 0-1200 (-180 to +180 degrees) in 0.3 degree steps",
            }
        )

    if "HorizontalRange" in ds:
        ds["HorizontalRange"].attrs.update(
            {
                "units": "m",
                "long_name": "horizontal distance along seabed from sonar to measurement point",
            }
        )

    if "OST_62" in ds:
        ds["OST_62"].attrs.update(
            {
                "units": "percent",
                "long_name": "Oxygen percent saturation",
                "standard_name": "fractional_saturation_of_oxygen_in_sea_water",
            }
        )

    if "P_1" in ds:
        ds["P_1"].attrs.update(
            {
                "units": "dbar",
                "long_name": "Uncorrected pressure",
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
            ds["P_1ac"].attrs.update({"note": ds.attrs["P_1ac_note"]})

    if "PAR_905" in ds:
        ds["PAR_905"].attrs.update(
            {
                "units": "umol m-2 s-1",
                "long_name": "Photosynthetically active " "radiation",
            }
        )

    if "pH_159" in ds:
        ds["pH_159"].attrs.update(
            {
                "units": "1",
                "standard_name": "sea_water_ph_reported_on_total_scale",
            }
        )

    if "Ptch_1216" in ds:
        ds["Ptch_1216"].attrs.update(
            {
                "units": "degrees",
                "long_name": "Instrument Pitch",
                "standard_name": "platform_pitch",
            }
        )

    # MIGHT NEED TO SPECIFY SONAR INSTRUMENT TYPE
    if "points" in ds:
        ds["points"].attrs.update(
            {
                "units": "1",
                "long_name": "point number in scan",
            }
        )

    if "ProfileRange" in ds:
        ds["ProfileRange"].attrs.update(
            {
                "units": "sample unit",
                "long_name": "First digitized range value above threshold in sample units",
                "comment": "Sample units are based on a sound velocity of 1500 m/s. For ranges <5m, one sample unit = 2 mm. For ranges >=5 m, one sampe unit = 10 mm.",
            }
        )

    if "RH_910" in ds:
        ds["RH_910"].attrs.update(
            {
                "units": "percent",
                "standard_name": "relative_humidity",
            }
        )

    if "Roll_1217" in ds:
        ds["Roll_1217"].attrs.update(
            {
                "units": "degrees",
                "long_name": "Instrument Roll",
                "standard_name": "platform_roll",
            }
        )

    if "Rn_963" in ds:
        ds["Rn_963"].attrs.update(
            {
                "units": "mm",
                "standard_name": "thickness_of_rainfall_amount",
            }
        )

    if "S_41" in ds:
        ds["S_41"].attrs.update(
            {
                "units": "1",
                "long_name": "Salinity, PSU",
                "comment": "Practical salinity units (PSU)",
                "standard_name": "sea_water_practical_salinity",
            }
        )

    if "sample" in ds:
        ds["sample"].attrs.update(
            {
                "units": "1",
                "long_name": "sample number",
            }
        )

    if "scan" in ds:
        ds["scan"].attrs.update(
            {
                "units": "1",
                "long_name": "scan number in sweep",
            }
        )

    if "SlantRange" in ds:
        ds["SlantRange"].attrs.update(
            {
                "units": "m",
                "long_name": "slant distance from sonar to seabed",
            }
        )

    if "sonar_hgt" in ds:
        ds["sonar_hgt"].attrs.update(
            {
                "units": "m",
                "long_name": "sonar height above seabed",
                "standard_name": "height_above_sea_floor",
            }
        )

    if "sonar_image" in ds:
        ds["sonar_image"].attrs.update(
            {
                "units": "1",
                "long_name": "sonar image data",
                "comment": "Units are integer values from 0-127",
                "standard_name": "acoustic_area_backscattering_strength_in_sea_water",
            }
        )

        if "sonar_image_note" in ds.attrs:
            utils.insert_note(ds, "sonar_image", ds.attrs["sonar_image_note"])

    if "SonarAngle" in ds:
        ds["SonarAngle"].attrs.update(
            {
                "units": "degrees",
                "long_name": "Instrument angle",
                "comment": "Angle = 0.3 x (Sonar Position - 600)",
            }
        )

    if "SonarPosition" in ds:
        ds["SonarPosition"].attrs.update(
            {
                "units": "1",
                "long_name": "Instrument position",
                "comment": "Integer values 0-1200 (-180 to +180 degrees) in 0.3 degree steps",
            }
        )

    if "SpC_48" in ds:
        ds["SpC_48"].attrs.update(
            {
                "units": "S/m",
                "long_name": "Specific Conductivity",
                "comment": "Temperature compensated to 25 °C",
                "standard_name": "sea_water_electrical_conductivity_at_reference_temperature",
            }
        )

    if "StepDirection" in ds:
        ds["StepDirection"].attrs.update(
            {
                "units": "1",
                "long_name": "Transducer head step direction",
                "comment": "0 = counter-clockwise, 1 = clockwise",
            }
        )

    if "sweep" in ds:
        ds["sweep"].attrs.update(
            {
                "units": "1",
                "long_name": "sweep number",
            }
        )

    if "T_21" in ds:
        ds["T_21"].attrs.update(
            {"units": "degree_C", "standard_name": "air_temperature"}
        )

    if "T_28" in ds:
        ds["T_28"].attrs.update(
            {
                "units": "degree_C",
                "standard_name": "sea_water_temperature",
                "long_name": "Temperature",
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

    if "TALPErfu" in ds:
        ds["TALPErfu"].attrs.update(
            {
                "units": "percent",
                "long_name": "Total algae phycoerythrin, RFU",
                "comment": "Relative fluorescence units (RFU); formerly called BGAPErfu (Blue green algae phycoerythrin, RFU)",
            }
        )

    if "theta" in ds:
        ds["theta"].attrs.update(
            {
                "units": "radians",
                "long_name": "Head angle relative to true north corrected for heading offset in north up convention",
                "comment": "Use theta and horizontal range to plot image data in polar convention",
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

    if "Turb_std" in ds:
        ds["Turb_std"].attrs.update(
            {
                "units": "1",
                "long_name": "Turbidity burst standard deviation (NTU)",
                "comment": "Nephelometric turbidity units (NTU)",
                "cell_methods": "time: standard_deviation",
            }
        )

    if "Tx_1211" in ds:
        ds["Tx_1211"].attrs.update(
            {
                "units": "degree_C",
                "long_name": "Instrument Internal Temperature",
            }
        )

    if "u_1205" in ds:
        ds["u_1205"].attrs.update(
            {
                "units": "m s^-1",
                "long_name": "Eastward Velocity",
            }
        )

    if "v_1206" in ds:
        ds["v_1206"].attrs.update(
            {
                "units": "m s^-1",
                "long_name": "Northward Velocity",
            }
        )

    if "water_level" in ds:
        ds["water_level"].attrs.update(
            {
                "units": "m",
                "long_name": "Water level NAVD88",
                "standard_name": "sea_surface_height_above_geopotential_datum",
                "geopotential_datum_name": "NAVD88",
            }
        )

    if "WD_410" in ds:
        ds["WD_410"].attrs.update(
            {
                "units": "degrees",
                "long_name": "mean wind from direction relative to true north",
                "standard_name": "wind_from_direction",
            }
        )

    if "WD_gust" in ds:
        ds["WD_gust"].attrs.update(
            {
                "units": "degrees",
                "long_name": "maximum wind from direction relative to true north",
                "standard_name": "wind_gust_from_direction",
            }
        )

    if "WD_min" in ds:
        ds["WD_min"].attrs.update(
            {
                "units": "degrees",
                "long_name": "minimum wind from direction relative to true north",
            }
        )

    if "WG_402" in ds:
        ds["WG_402"].attrs.update(
            {
                "units": "m/s",
                "long_name": "maximum wind speed",
                "standard_name": "wind_speed_of_gust",
            }
        )

    if "WS_401" in ds:
        ds["WS_401"].attrs.update(
            {
                "units": "m/s",
                "long_name": "mean wind speed",
                "standard_name": "wind_speed",
            }
        )

    if "WS_min" in ds:
        ds["WS_min"].attrs.update({"units": "m/s", "long_name": "minimum wind speed"})

    if "x" in ds:
        ds["x"].attrs.update(
            {
                "units": "m",
                "long_name": "easterly horizontal distance from sonar",
            }
        )

    if "y" in ds:
        ds["y"].attrs.update(
            {
                "units": "m",
                "long_name": "northerly horizontal distance from sonar",
            }
        )

    return ds
