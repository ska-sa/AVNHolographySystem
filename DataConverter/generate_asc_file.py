import h5py
import numpy as np
import pandas as pd
import sys
import katpoint
import scipy.interpolate as interpolate

# Filenames from arguments
h5_filename = sys.argv[1]
csv_filename = sys.argv[2]
snp_filename = sys.argv[3]

# Let's make the pmodl file optional...
try:
    pmodl_filename = sys.argv[4]
    pmodl_file = open(pmodl_filename, "r")
except IndexError:
    print "Warning! No pmodl file supplied, using a default zero one."
    pmodl_file = None
except IOError:
    print "Error! %s doesn't exist, using zero pointing model." % pmodl_filename
    pmodl_file = None
h5_file = h5py.File(h5_filename)
csv_file = pd.read_csv(csv_filename, skipinitialspace=True)

#TODO: make this not a raw_input thing, it's messy.
channel_num = int(raw_input("Please enter the channel number: "))
freq = float(raw_input("Enter the frequency (MHz) of the target: "))
pol_selection = int(raw_input("1 for AUTxRef, 2 for RefXAUT: "))

data_column = np.array(h5_file["Data/VisData"][:,channel_num,pol_selection,0]
                       + 1j*h5_file["Data/VisData"][:,channel_num,pol_selection,1])

csv_timestamps = np.array(csv_file["Timestamp"])
h5_timestamps = np.array(h5_file["Data/Timestamps"])


# This cut straight from my previous Katdal work.
# Should probably get around to modularising that properly.
# Credit to http://stackoverflow.com/a/12737895 for this function:
def decdeg2dms(dd):
    negative = dd < 0
    dd = abs(dd)
    minutes, seconds = divmod(dd * 3600, 60)
    degrees, minutes = divmod(minutes, 60)
    if negative:
        if degrees > 0:
            degrees = -degrees
        elif minutes > 0:
            minutes = -minutes
        else:
            seconds = -seconds
    return int(degrees), int(minutes), float(seconds)


def fs_to_kp_pointing_model(pmodl_file):
    """Parses a Field System pointing model file (mdlpo.ctl)
    Returns:
        katpoint pointing model (object? string?)
    """
    if pmodl_file is None:
        pmodl_string = "0 " * 22
        pmodl_string = pmodl_string[:-1]  # Remove the resulting space on the end.
        return pmodl_string

    if not isinstance(pmodl_file, file):
        raise TypeError("%s not a text file." % (repr(pmodl_file)))
    lines = []
    for i in pmodl_file:
        lines.append(i)
    del i
    if len(lines) != 19:
        raise TypeError(
            "%s not correct length for pointing model file. File is %d lines long." % (repr(pmodl_file), len(lines)))
    # Line 5 gives the enabled parameters:
    params_implemented = lines[5].split()

    if len(params_implemented) != 31:
        raise TypeError("%s not correct format for pointing model file." % (repr(pmodl_file)))
    # The first number on the line is the phi value
    phi = float(params_implemented[0])

    # If any of the higher ones are used, throw a warning:
    # TODO: figure out how to do this properly in Python. There must be a prettier way.
    if params_implemented[23] == '1' or \
            params_implemented[24] == '1' or \
            params_implemented[25] == '1' or \
            params_implemented[26] == '1' or \
            params_implemented[27] == '1' or \
            params_implemented[28] == '1' or \
            params_implemented[29] == '1' or \
            params_implemented[30] == '1':
        print "Warning: params 23 - 30 are not used in katpoint, but are used in the pointing model file."

    # Lines 7, 9, 11, 13, 15 and 17 each have 5 parameters on them.
    params = [0]  # Pad place number 0 so that the numbers correspond.
    params.extend(lines[7].split())
    params.extend(lines[9].split())
    params.extend(lines[11].split())
    params.extend(lines[13].split())
    params.extend(lines[15].split())
    params.extend(lines[17].split())

    for i in range(len(params)):
        params[i] = float(params[i])

    pmodl_string = ""

    for i in range(1, 23):
        if params_implemented[i] == '1' and float(params[i]) != 0:
            pmodl_string += "%02d:%02d:%06.3f " % (decdeg2dms(
                float(params[i])))  # I had thought that Katpoint needs dd:mm:ss.xx format, but apparently it doesn't.
            # pmodl_string += "%s "%(params[i])

        else:
            pmodl_string += "0 "
    del i

    pmodl_string = pmodl_string[:-1]  # Remove the resulting space on the end.
    return pmodl_string, np.array(params[1:])


pmodl_set = fs_to_kp_pointing_model(pmodl_file)
# TODO: Update this for whichever AVN antenna is being used.
antenna_str = "ant1, 5:45:2.48, -0:18:17.92, 116, 32.0, 0 0 0, %s" % (pmodl_set[0])
pmodl_set = np.array(pmodl_set[1:]).transpose()
antenna = katpoint.Antenna(antenna_str)
satellite_name = "Eutelsat 7A"
TLE = """EUTELSAT 7A
1 28187U 04008A   18051.18726536  .00000042  00000-0  00000+0 0  9996
2 28187   0.0701 342.3392 0003208 340.3981 261.7175  1.00271618 51105"""
target = katpoint.Target("%s, tle, %s" % (satellite_name, TLE))
target.antenna = antenna

# TODO: check for limits.

# Interpolate CSV data onto RF data:
az_interpolator = interpolate.interp1d(csv_timestamps, np.array(csv_file["Azim actual position"]), kind="cubic")
csv_az = az_interpolator(h5_timestamps)
el_interpolator = interpolate.interp1d(csv_timestamps, np.array(csv_file["Elev actual position"]), kind="cubic")
csv_el = el_interpolator(h5_timestamps)

# TODO: Flagging. Don't have any of that yet.

# Find the right lines in the csv log file.
with open("output_file.asc", "w") as output_file:
    output_file.write("#freq_MHz=%.4f\n" % freq)
    output_file.write("#targetaz_deg=%.4f\n" % target_az)
    output_file.write("#targetel_deg=%.4f\n" % target_el)
    for i in range(len(h5_timestamps)):
        azel = np.degrees(antenna.pointing_model.reverse(np.radians(csv_az[i]),
                                                         np.radians(csv_el[i])))
        output_file.write("%.6f\t%.6f\t%.6f\t%.6f\t%.2f\n" %
                          (azel[0] - target_az, azel[1] - target_el,
                           np.abs(data_column[i]), np.degrees(np.angle(data_column[i])), h5_timestamps[i]))
