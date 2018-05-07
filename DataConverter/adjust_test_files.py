import h5py
import numpy as np
import time
import pandas as pd

# copied from snp file. minus and plus later are to ensure that the intervals are longer.
snp_start_time = 1519912800
snp_stop_time = 1519963830

adjust_h5 = False
adjust_csv = True

if adjust_h5:
    h5_start_time = snp_start_time + 10
    h5_stop_time = snp_stop_time - 10

    h5_file = h5py.File("test.h5", "r+")
    h5_time_length = len(h5_file["Data/Timestamps"])

    new_h5_timestamps = np.linspace(h5_start_time, h5_stop_time, h5_time_length, endpoint=True, dtype=np.float64)
    del h5_file["Data/Timestamps"]
    h5_file["Data/Timestamps"] = new_h5_timestamps

    h5_file.close()

if adjust_csv:
    csv_start_time = snp_start_time - 10
    csv_stop_time = snp_stop_time + 10

    csv_file = pd.read_csv("test.csv", skipinitialspace=True)

    csv_time_length = len(csv_file["Timestamp"])
    new_csv_timestamps = np.linspace(csv_start_time, csv_stop_time, csv_time_length, dtype=np.float64)
    new_csv_timestamps *= 1e3
    new_csv_timestamps = new_csv_timestamps.astype(np.int)
    csv_file["Timestamp"] = new_csv_timestamps

    csv_file.to_csv("test_out.csv", index=False)