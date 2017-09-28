import h5py
import numpy as np
import pandas as pd
import sys

h5_filename = sys.argv[1]
csv_filename = sys.argv[2]
h5_file = h5py.File(h5_filename)
csv_file = pd.read_csv(csv_filename, skipinitialspace=True)

channel_num = int(raw_input("Please enter the channel number: "))
pol_selection = int(raw_input("1 for AUTxRef, 2 for RefXAUT: "))

data_column = np.array(h5_file["Data/VisData"][:,channel_num,pol_selection,0] + 1j*h5_file["Data/VisData"][:,channel_num,pol_selection,1])

csv_timestamps = np.array(csv_file["Timestamp"])
h5_timestamps = np.array(h5_file["Data/Timestamps"])

# Find the right lines in the csv log file.
with open("output_file.asc", "w") as output_file:
    output_file.write()
    for i in range(len(h5_timestamps)):
        idx = np.abs(csv_timestamps - h5_timestamps[i]).argmin()
        output_file.write("%.6f\t%.6f\t%.6f\t%.6f\t%.2f\n" %
                          (csv_file["Azim actual position"][idx], csv_file["Elev actual position"][idx],
                           np.abs(data_column[i]), np.degrees(np.angle(data_column[i])), h5_timestamps[i]))
