# TODO update this for Jon's Q&D log file.

import h5py
import numpy as np
import pandas as pd
import sys
import scipy.interpolate as interpolate
from astropy.time import Time
import time


h5_filename = sys.argv[1]
csv_filename = sys.argv[2]

h5_file = h5py.File(h5_filename)
csv_file = pd.read_csv(csv_filename, skipinitialspace=True)

channel_num = int(raw_input("Please enter the channel number: "))
freq = float(raw_input("Enter the frequency (MHz) of the target: "))
pol_selection = int(raw_input("1 for AUTxRef, 2 for RefXAUT: "))
target_ha = float(raw_input("Enter target hour angle (degrees): "))
target_dec = float(raw_input("Enter target declination (degrees): "))

data_column = np.array(h5_file["Data/VisData"][:,channel_num,pol_selection,0] + 1j*h5_file["Data/VisData"][:,channel_num,pol_selection,1])

# Convert from MJD to Unix time
csv_mjd = np.array(csv_file["# MJD"])
csv_timestamps = np.zeros(len(csv_mjd), dtype=np.float)
for i in range(len(csv_mjd)):
    t = Time(csv_mjd[i], format="mjd")
    csv_timestamps[i] = t.unix
print "CSV start time: %s" % (time.strftime("%Y.%m.%D %H:%M:%S",time.gmtime(csv_timestamps[0])))
print "CSV end time: %s" % (time.strftime("%Y.%m.%D %H:%M:%S",time.gmtime(csv_timestamps[-1])))
# Timestamps recorded in the h5 at the *start* of the accumulation, I think, so I'm going to
# shift them forward by 1/2 an accumulation length in order to get them in the average
# in order to better line up the pointing.
h5_timestamps = np.array(h5_file["Data/Timestamps"])
h5_timestamps += (h5_timestamps[1] - h5_timestamps[0])/2
print "H5 start time: %s" % (time.strftime("%Y.%m.%D %H:%M:%S",time.gmtime(h5_timestamps[0])))
print "H5 end time: %s" % (time.strftime("%Y.%m.%D %H:%M:%S",time.gmtime(h5_timestamps[-1])))

# Timestamps won't line up exactly so we'll interpolate the angle information.
hourangle_interpolator = interpolate.interp1d(csv_timestamps, np.array(csv_file["HA"]), kind="cubic")
new_hourangle = hourangle_interpolator(h5_timestamps)
declination_interpolator = interpolate.interp1d(csv_timestamps, np.array(csv_file["Dec"]), kind="cubic")
new_declination = declination_interpolator(h5_timestamps)

# Find the right lines in the csv log file.
with open("output_file.asc", "w") as output_file:
    output_file.write("#freq_MHz=%.4f\n" % freq)
    output_file.write("#targetha_deg=%.4f\n" % target_ha)
    output_file.write("#targetdec_deg=%.4f\n" % target_dec)
    for i in range(len(h5_timestamps)):
        output_file.write("%.6f\t%.6f\t%.6f\t%.6f\t%.2f\n" %
                          (new_hourangle[i] - target_ha, new_declination[i] - target_dec,
                           np.abs(data_column[i]), np.degrees(np.angle(data_column[i])), h5_timestamps[i]))
