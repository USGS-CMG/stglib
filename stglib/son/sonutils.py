import datetime
import struct

from ..core import utils


def bin2ascii(binString):
    """
    Read a set of binary ASCII values and return a string
    Input:
        binString: Binary ASCII value(s)
    Output:
        output: String of characters from binString
    """
    # Create an array of characters
    charArray = [chr(x) for x in binString]

    # Convert the array of characters into one string
    output = ""
    for char in charArray:
        output += char

    # Remove ASCII null values (00 in hex) and return the string
    return output.replace("\x00", "")


def parse_pingHeader(fheader, month_day):
    """
    Read the universal 81R header using the guide provided by the manual
    Inputs:
        fheader: The file header in binary form
        timeprefix: a four digit date string in the format of mmdd extracted
            from the filename
    Outputs:
        time: The time at the current ping
        pingheader: A dictionary containing information extracted from the header
    """

    # First 3 bytes read 81R
    pingheader = {}
    pingheader["SONReturnDataHeaderType"] = bin2ascii(fheader[0:3])

    # Get sonar model
    models = ["881L-GS", "881A-GS", "882L", "882A"]
    son_type = fheader[3]
    pingheader["SONSonarType"] = "Imagenex " + models[son_type]

    # Get file size information based on the model of the sonar
    switch_size = [128, 40]
    pingheader["SONSwitchCommandBytes"] = switch_size[fheader[3] % 2]
    return_size = [256, 32]
    pingheader["SONReturnHeaderBytes"] = return_size[fheader[3] % 2]
    pingheader["SONTotalBytes"] = struct.unpack("<I", fheader[4:8])[0]

    # Extract information
    tx_orient = ["Down", "Up"]
    orientation = int(bin(fheader[319])[2:].zfill(8)[7])  # 1 is up
    pingheader["SONOrientation"] = tx_orient[orientation]
    modes = ["SONSector", "Polar", "Sidescan"]
    pingheader["SONMode"] = modes[fheader[324]]
    pingheader["SONRangeOffset"] = struct.unpack("<f", fheader[325:329])[0]
    pingheader["SONSoundVelocity"] = struct.unpack("<f", fheader[338:342])[0]
    pingheader["SONTransmitFrequency"] = struct.unpack("<f", fheader[342:346])[0]
    pingheader["SONPingRepetitionRate"] = struct.unpack("<f", fheader[346:350])[0]
    pingheader["SONSamplesPerPing"] = struct.unpack("<L", fheader[353:357])[0]
    # pingheader["SONAcousticRangeSetting"] = struct.unpack("<f", fheader[369:373])[0]
    pingheader["SONRangeResolution"] = struct.unpack("<f", fheader[373:377])[0]
    pingheader["SONPingNumber"] = struct.unpack("<L", fheader[377:381])[0]
    pingheader["SONSystemFlag"] = struct.unpack("<B", fheader[381:382])[0]
    pingheader["SONGyroStatus"] = struct.unpack("<B", fheader[382:383])[0]
    pingheader["SONMountingAngleOffset"] = struct.unpack("<f", fheader[383:387])[0]
    pingheader["SONLocalLatitude"] = struct.unpack("<f", fheader[387:391])[0]
    pingheader["SONCompassDeclination"] = struct.unpack("<f", fheader[391:395])[0]

    # Concatenate date
    dstr = month_day + bin2ascii(fheader[14:27])
    in_fmt = "%m%d%Y%H%M%S.%f"
    time = datetime.datetime.strptime(dstr, in_fmt)
    return time, pingheader


def parse_switchCommand(scommand):
    """
    Parse the switch commands unique to the 881A-GS Model
    Input:
        scommand: Switch commands in binary format
    Output:
        SwitchCommand: A dictionary containing information extracted
    """

    # All conversions are found in the manual
    SwitchCommand = {}
    SwitchCommand["SONHeadID"] = scommand[2]
    SwitchCommand["SONRange"] = scommand[3]
    SwitchCommand["SONStartGain"] = scommand[8]
    SwitchCommand["SONLogF"] = 10 * (scommand[9] + 1)
    SwitchCommand["SONAbsorption"] = scommand[10] / 100
    SwitchCommand["SONTrainAngle"] = 3 * scommand[11] - 180
    SwitchCommand["SONSectorWidth"] = 3 * scommand[12]
    SwitchCommand["SONStepSize"] = 0.3 * scommand[13]
    SwitchCommand["SONPulseLength"] = 10 * scommand[14]
    SwitchCommand["SONMinRange"] = scommand[15] / 10
    SwitchCommand["SONNDataPoints"] = scommand[19] * 10
    SwitchCommand["SONDataBits"] = scommand[20]
    profile_cmd = ["OFF", "ON"]
    profile = scommand[23]
    SwitchCommand["SONProfile"] = profile_cmd[profile]
    SwitchCommand["SONFrequency"] = 175 + scommand[25] * 5

    return SwitchCommand


def parse_returnHeader(sheader, switches):
    """
    Parse the return data header for the 881A-GS model sonar
    Inputs:
        sheader: The return data header in binary form
        switches: A dictionary of switch commands created using parseSwitchCommand()
    Output:
        returnHeader: A dictionary containing data extracted from the return data header
    """

    # A 3-character string that used to determine the size of the sonar output
    returnHeader = {}
    returnHeader["ReturnDataHeaderType"] = bin2ascii(sheader[0:3])

    # ID of the transducer head
    # returnHeader["HeadID"] = hex(sheader[3])

    # Extracted using Doug Wilson's method
    returnHeader["HeadPosition"] = (63 & sheader[6]) * 128 + (127 & sheader[5])
    returnHeader["HeadAngle"] = (returnHeader["HeadPosition"] - 600) * switches[
        "SONStepSize"
    ]

    # Extract additional information
    # Method for step direction taken from the manual
    returnHeader["StepDirection"] = (sheader[6] & 64) >> 6
    # returnHeader["Range"] = sheader[7]

    # The following variables are all extracted using methods provided by the manual
    # Profile Range
    HB = (sheader[9] & 0x7E) >> 1
    LB = ((sheader[9] & 0x01) << 7) | (sheader[8] & 0x7F)
    returnHeader["ProfileRange"] = (HB << 8) | LB

    # Data Bytes
    HB = (sheader[11] & 0x7E) >> 1
    LB = ((sheader[11] & 0x01) << 7) | (sheader[10] & 0x7F)
    returnHeader["NDataBytes"] = (HB << 8) | LB

    # Sonar Position
    HB = (sheader[13] & 0x7E) >> 1
    LB = ((sheader[13] & 0x01) << 7) | (sheader[12] & 0x7F)
    returnHeader["SonarPosition"] = (HB << 8) | LB
    returnHeader["SonarAngle"] = 0.3 * (returnHeader["SonarPosition"] - 600)

    # Pitch
    HB = (sheader[15] & 0x7E) >> 1
    LB = ((sheader[15] & 0x01) << 7) | (sheader[14] & 0x7F)
    returnHeader["Pitch"] = (
        (((HB << 8) | LB) - 16384 * (int(bin(sheader[15])[-1]))) * 360 / 16384
    )

    # Roll
    HB = (sheader[17] & 0x7E) >> 1
    LB = ((sheader[17] & 0x01) << 7) | (sheader[16] & 0x7F)
    returnHeader["Roll"] = ((HB << 8) | LB) * 360 / 16384

    # Heading
    HB = (sheader[19] & 0x7E) >> 1
    LB = ((sheader[19] & 0x01) << 7) | (sheader[18] & 0x7F)
    returnHeader["Heading"] = ((HB << 8) | LB) * 360 / 16384

    # Gyro Heading
    HB = (sheader[22] & 0x7E) >> 1
    LB = ((sheader[22] & 0x01) << 7) | (sheader[21] & 0x7F)
    returnHeader["GyroHeading"] = ((HB << 8) | LB) * 360 / 16384

    # Use the previously extracted 3-character string and data from SwitchCommand to determine
    # the amount of data points per ping. Table containing this information is in the manual
    if returnHeader["ReturnDataHeaderType"] == "INA":
        if switches["SONDataBits"] == 4:
            returnHeader["NReturnBytes"] = 128
            returnHeader["NPoints"] = 256
        elif switches["SONDataBits"] == 8:
            returnHeader["NReturnBytes"] = 252
            returnHeader["NPoints"] = 252
        elif switches["SONDataBits"] == 16:
            returnHeader["NReturnBytes"] = 500
            returnHeader["NPoints"] = 250
    elif returnHeader["ReturnDataHeaderType"] == "INB":
        if switches["SONDataBits"] == 4:
            returnHeader["NReturnBytes"] = 252
            returnHeader["NPoints"] = 504
        elif switches["SONDataBits"] == 8:
            returnHeader["NReturnBytes"] = 500
            returnHeader["NPoints"] = 500
        elif switches["SONDataBits"] == 16:
            returnHeader["NReturnBytes"] = 500
            returnHeader["NPoints"] = 250
    elif returnHeader["ReturnDataHeaderType"] == "INC":
        returnHeader["NPoints"] = 0
    else:
        print(
            f"{returnHeader['ReturnDataHeaderType']} is not a recognized data header type"
        )

    return returnHeader


def ds_add_attrs(ds):
    """
    Add attributes: units, standard name from CF website, long names
    """

    ds["time"].attrs.update(
        {"standard_name": "time", "axis": "T", "long_name": "time (UTC)"}
    )

    if "sweep" in ds:
        ds["sweep"].attrs.update(
            {
                "units": "1",
                "long_name": "sweep number",
            }
        )

    if "scan" in ds:
        ds["scan"].attrs.update(
            {
                "units": "1",
                "long_name": "scan number in sweep",
            }
        )

    if "points" in ds:
        ds["points"].attrs.update(
            {
                "units": "1",
                "long_name": "point number in scan",
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

    if "Roll_1217" in ds:
        ds["Roll_1217"].attrs.update(
            {
                "units": "degrees",
                "long_name": "Instrument Roll",
                "standard_name": "platform_roll",
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

    if "GyroHeading" in ds:
        ds["GyroHeading"].attrs.update(
            {
                "units": "degrees",
                "long_name": "Gyro Heading",
                "standard_name": "platform_orientation",
            }
        )

    if "sonar_image" in ds:
        ds["sonar_image"].attrs.update(
            {
                "units": "1",
                "long_name": "sonar image data",
                "comments": "Units are integer values from 0-127",
                "standard_name": "acoustic_area_backscattering_strength_in_sea_water",
            }
        )

        if "sonar_image_note" in ds.attrs:
            utils.insert_note(ds, "sonar_image", ds.attrs["sonar_image_note"])

    if "sonar_hgt" in ds:
        ds["sonar_hgt"].attrs.update(
            {
                "units": "m",
                "long_name": "sonar height above seabed",
                "standard_name": "height_above_sea_floor",
            }
        )

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

    if "HeadPosition" in ds:
        ds["HeadPosition"].attrs.update(
            {
                "units": "1",
                "long_name": "Transducer head position",
                "comments": "Integer values 0-1200 (-180 to +180 degrees) in 0.3 degree steps",
            }
        )

    if "HeadAngle" in ds:
        ds["HeadAngle"].attrs.update(
            {
                "units": "degrees",
                "long_name": "Transducer head angle",
                "comments": "Angle = 0.3 x (Head Position - 600)",
            }
        )

    if "StepDirection" in ds:
        ds["StepDirection"].attrs.update(
            {
                "units": "1",
                "long_name": "Transducer head step direction",
                "comments": "0 = counter-clockwise, 1 = clockwise",
            }
        )

    if "ProfileRange" in ds:
        ds["ProfileRange"].attrs.update(
            {
                "units": "sample unit",
                "long_name": "First digitized range value above threshold in sample units",
                "comments": "Sample units are based on a sound velocity of 1500 m/s. For ranges <5m, one sample unit = 2 mm. For ranges >=5 m, one sampe unit = 10 mm.",
            }
        )

    if "SonarPosition" in ds:
        ds["SonarPosition"].attrs.update(
            {
                "units": "1",
                "long_name": "Instrument position",
                "comments": "Integer values 0-1200 (-180 to +180 degrees) in 0.3 degree steps",
            }
        )

    if "SonarAngle" in ds:
        ds["SonarAngle"].attrs.update(
            {
                "units": "degrees",
                "long_name": "Instrument angle",
                "comments": "Angle = 0.3 x (Sonar Position - 600)",
            }
        )

    if "SlantRange" in ds:
        ds["SlantRange"].attrs.update(
            {
                "units": "m",
                "long_name": "slant distance from sonar to seabed",
            }
        )

    if "HorizontalRange" in ds:
        ds["HorizontalRange"].attrs.update(
            {
                "units": "m",
                "long_name": "horizontal distance along seabed from sonar to measurement point",
            }
        )

    if "theta" in ds:
        ds["theta"].attrs.update(
            {
                "units": "radians",
                "long_name": "Head angle relative to true north corrected for heading offset in north up convention",
                "comments": "Use theta and horizontal range to plot image data in polar convention",
            }
        )

    return ds
