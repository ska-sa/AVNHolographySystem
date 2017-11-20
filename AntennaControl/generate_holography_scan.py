import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import time
import katpoint
from optparse import OptionParser



def generate_raster(total_extent_az, total_extent_el, az_resolution=0.05, el_resolution=0.05, dwell_time=10,
                    n_slew_points=5, slew_speed=0.02, scan_speed=0.01):
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

    for i in range(n_scans):
        # Time: the division is making a lot of zeros.
        # Slew to the start of the next scan
        x_points = np.concatenate((x_points, np.linspace(x_points[-1], x_space[0], n_slew_points + 1, endpoint=False)[1:]))
        y_points = np.concatenate((y_points, np.linspace(y_points[-1], y_space[i], n_slew_points + 1, endpoint=False)[1:]))
        time_points = np.concatenate((time_points, np.ones(n_slew_points, dtype=int)*float(slew_time)/n_slew_points))
        for j in range(n_slew_points):
            label.append("slew")

        # Scan
        x_points = np.concatenate((x_points, x_space))
        y_points = np.concatenate((y_points, np.ones(len(x_space))*y_space[i]))
        time_points = np.concatenate((time_points, np.ones(len(x_space), dtype=int) * float(scan_time) / n_points_per_scan))
        for j in range(len(x_space)):
            label.append("scan")

        # Slew back to the target
        x_points = np.concatenate((x_points, np.linspace(x_points[-1], 0, n_slew_points + 1, endpoint=False)[1:]))
        y_points = np.concatenate((y_points, np.linspace(y_points[-1], 0, n_slew_points + 1, endpoint=False)[1:]))
        time_points = np.concatenate((time_points, np.ones(n_slew_points, dtype=int) * float(slew_time) / n_slew_points))
        for j in range(n_slew_points):
            label.append("slew")

        x_points = np.concatenate((x_points, [0]))
        y_points = np.concatenate((y_points, [0]))
        time_points = np.concatenate((time_points, [dwell_time]))
        label.append("track")

    time_points = np.cumsum(time_points)
    time_points -= time_points[0]

    return x_points, y_points, time_points, label


def update_line(num, data, line):
    line.set_data(data[..., :num])
    return line,


if __name__ == "__main__":
    parser = OptionParser(usage="usage: %prog <AZ> <EL> <START_DATE> <START_TIME> [options]", version="%prog 0.1")
    parser.add_option("-x", "--x-size", action="store", type=float, dest="raster_size_az", default=1.0,
                      help="Extent (degrees) of raster scan in azimuth. Default: 1.0")
    parser.add_option("-y", "--y-size", action="store", type=float, dest="raster_size_el", default=1.0,
                      help="Extent (degrees) of raster scan in elevation. Default: 1.0")
    parser.add_option("-o", "--output-filename", action="store", type=str, dest="output_filename", default="output",
                      help="Name of the snap file to be generated. Default: output")
    parser.add_option("-p", "--plot", action="store_true", dest="plot", default=False,
                      help="Plot animated illustration of scan on an x-y grid.")

    options, args = parser.parse_args()

    if len(args) != 4:
        parser.error("Incorrect number of arguments! Expected 4: <AZ> <EL> <START_DATE> <START_TIME>")
    target_az, target_el = float(args[0]), float(args[1])

    raster_size_az = options.raster_size_az
    raster_size_el = options.raster_size_el

    output_filename = options.output_filename
    start_time_input = args[2] + " " + args[3]
    try:
        start_time = time.mktime(time.strptime(start_time_input, "%Y-%m-%d %H:%M:%S"))
    except ValueError:
        print "Either start date or time are incorrectly formatted. Use format: 2001-09-11 09:57:55"
        exit(-1)

    myTarget = katpoint.construct_azel_target(np.radians(target_az), np.radians(target_el))
    antenna_str = "Kuntunse, 5:45:2.48, -0:18:17.92, 116, 32.0"
    antenna = katpoint.Antenna(antenna_str)
    myTarget.antenna = antenna
    d_az, d_el, t, comment = generate_raster(raster_size_az, raster_size_el)

    #print len(d_az), len(t), len(comment)

    with open("%s.snp" % output_filename, "w") as output_file:
        output_file.write("\" Holography raster script developed 2017-07-19 by James Smith.\n")
        output_file.write("\" Field system snap file automagically generated at: %s." % time.strftime("%Y-%m-%d %H:%M:%S"))
        output_file.write("\" File: %s.snp\n" % output_filename)
        output_file.write("\" Target Az/El: %.2f, %.2f\n" % (target_az, target_el))
        output_file.write("\" Raster size: %.2f x %.2f degrees\n" % (raster_size_az, raster_size_el))
        output_file.write("\" Schedule to start at: %s\n" % start_time_input)  # TODO: This needs fixing.
        for i in range(len(d_az)):
            output_file.write("\"#%d,%s,nominal\n!%s\nsource=azel,%.6fd,%.6fd\n" %
                              (start_time + t[i], comment[i],
                               time.strftime("%Y.%j.%H:%M:%S",time.gmtime(start_time + t[i])),
                               target_az + d_az[i], target_el + d_el[i]))
        output_file.write("stop")

    if options.plot:
        x, y = np.degrees(myTarget.sphere_to_plane(np.radians(target_az + d_az), np.radians(target_el + d_el)))
        data = np.stack((x, y))
        fig = plt.figure()
        plt.xlim(np.min(x) - 0.5, np.max(x) + 0.5)
        plt.ylim(np.min(x) - 0.5, np.max(x) + 0.5)
        my_line, = plt.plot([], [], '.')
        line_ani = anim.FuncAnimation(fig, update_line, len(x), fargs=(data, my_line), interval=50, blit=True, repeat=False)

        line_ani.save("holog_scan.gif", writer="imagemagick")
        plt.show()
