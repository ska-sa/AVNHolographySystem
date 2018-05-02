import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import time
import calendar # for the one timegm function...
import katpoint

# Edit these before each run.
satellite_name = "Eutelsat 7A"
TLE = """EUTELSAT 7A
1 28187U 04008A   18051.18726536  .00000042  00000-0  00000+0 0  9996
2 28187   0.0701 342.3392 0003208 340.3981 261.7175  1.00271618 51105"""
antenna_str = "Kuntunse, 5:45:2.48, -0:18:17.92, 116, 32.0"
raster_size_az = 4  # degrees
raster_size_el = 4  # degrees
start_time_input = "2018-03-01 14:00:00"
output_filename = "output"


def generate_raster(total_extent_az, total_extent_el, az_resolution=0.05, el_resolution=0.05, dwell_time=10,
                    n_slew_points=5, slew_speed=0.02, scan_speed=0.01):
    """All angle values in degrees, speed values in degrees per second.
    """

    # These need to be forced as floats, otherwise integer division can bite you in the bum.
    total_extent_az = float(total_extent_az)
    total_extent_el = float(total_extent_el)
    slew_time = float(total_extent_az) / 2 / slew_speed
    scan_time = float(total_extent_az) / scan_speed
    n_scans = int(np.ceil(total_extent_el / el_resolution) + 1)
    n_points_per_scan = int(np.ceil(total_extent_az / az_resolution) + 1)
    if n_scans % 2 == 0:
        n_scans += 1

    y_space = np.linspace(-total_extent_el/2, total_extent_el/2, n_scans, endpoint=True)
    x_space = np.linspace(-total_extent_az/2, total_extent_az/2, n_points_per_scan, endpoint=False)

    x_points = np.array([0])
    y_points = np.array([0])
    time_points = np.array([dwell_time])
    label = ["track"]

    for ui in range(n_scans):
        # Time: the division is making a lot of zeros.
        # Slew to the start of the next scan
        x_points = np.concatenate((x_points,
                                   np.linspace(x_points[-1], x_space[0], n_slew_points + 1, endpoint=False)[1:]))
        y_points = np.concatenate((y_points,
                                   np.linspace(y_points[-1], y_space[ui], n_slew_points + 1, endpoint=False)[1:]))
        time_points = np.concatenate((time_points, np.ones(n_slew_points, dtype=int)*float(slew_time)/n_slew_points))
        for j in range(n_slew_points):
            label.append("slew")

        # Scan
        x_points = np.concatenate((x_points, x_space))
        y_points = np.concatenate((y_points, np.ones(len(x_space))*y_space[ui]))
        time_points = np.concatenate((time_points,
                                      np.ones(len(x_space), dtype=int) * float(scan_time) / n_points_per_scan))
        for j in range(len(x_space)):
            label.append("scan")

        # Slew back to the target
        x_points = np.concatenate((x_points, np.linspace(x_points[-1], 0, n_slew_points + 1, endpoint=False)[1:]))
        y_points = np.concatenate((y_points, np.linspace(y_points[-1], 0, n_slew_points + 1, endpoint=False)[1:]))
        time_points = np.concatenate((time_points,
                                      np.ones(n_slew_points, dtype=int) * float(slew_time) / n_slew_points))
        for j in range(n_slew_points):
            label.append("slew")

        x_points = np.concatenate((x_points, [0]))
        y_points = np.concatenate((y_points, [0]))
        time_points = np.concatenate((time_points, [dwell_time]))
        label.append("track")

    time_points = np.cumsum(time_points)
    time_points -= time_points[0]

    return x_points, y_points, time_points, label


def update_line(num, line_data, line):
    line.set_data(line_data[..., :num])
    return line,


if __name__ == "__main__":

    # target_az, target_el = float(args[0]), float(args[1])

    try:
        start_time = calendar.timegm(time.strptime(start_time_input, "%Y-%m-%d %H:%M:%S"))
    except ValueError:
        print "Either start date or time are incorrectly formatted. Use format: 2001-09-11 09:57:55"
        exit(-1)


    myTarget = katpoint.Target("%s, tle, %s" % (satellite_name, TLE))

    # noinspection PyTypeChecker
    # myTarget = katpoint.construct_azel_target(np.radians(target_az), np.radians(target_el))
    antenna = katpoint.Antenna(antenna_str)
    myTarget.antenna = antenna
    d_az, d_el, t, comment = generate_raster(raster_size_az, raster_size_el)

    x, y = [], []

    with open("%s.snp" % output_filename, "w") as output_file:
        output_file.write("\" Holography raster script developed 2018-02-21 by James Smith.\n")
        output_file.write("\" Field system snap file automagically generated at: %s."
                          % time.strftime("%Y-%m-%d %H:%M:%S"))
        output_file.write("\" File: %s.snp\n" % output_filename)
        output_file.write("\" Target: %s\n" % repr(myTarget))
        output_file.write("\" Target nominal az/el: %.2f, %.2f\n" % (np.degrees(myTarget.azel(start_time)[0]),
                                                                     np.degrees(myTarget.azel(start_time)[1])))
        output_file.write("\" Raster size: %.2f x %.2f degrees\n" % (raster_size_az, raster_size_el))
        output_file.write("\" Schedule to start at: %s\n" % start_time_input)
        for i in range(len(d_az)):
            output_file.write("\"#%d,%s,nominal\n!%s\nsource=azel,%.6fd,%.6fd\n" %
                              (start_time + t[i], comment[i],
                               time.strftime("%Y.%j.%H:%M:%S", time.gmtime(start_time + t[i])),
                               np.degrees(myTarget.azel(start_time + t[i])[0]) + d_az[i],
                               np.degrees(myTarget.azel(start_time + t[i])[1]) + d_el[i]))
            xy = myTarget.sphere_to_plane(np.radians(myTarget.azel(start_time + t[i])[0] + d_az[i]),
                                          np.radians(myTarget.azel(start_time + t[i])[1] + d_el[i]))
            x.append(np.degrees(xy[0]))
            y.append(np.degrees(xy[1]))
        output_file.write("stop")
