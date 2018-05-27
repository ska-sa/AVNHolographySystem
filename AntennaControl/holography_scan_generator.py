import numpy as np
import matplotlib.pyplot as plt
import time
import calendar  # for the one timegm function...
import katpoint

# TODO consider adding type checking, force the floats to be floats for example because otherwise the divisions will go wonky.
# TODO: Think about data validation. Scans-per-boresight should be less than overall scans, for example.

# Edit these before each run.
satellite_name = "SES-5"
TLE = """SES-5
1 38652U 12036A   18146.57916431  .00000021  00000-0  00000-0 0  9995
2 38652   0.0240 269.8242 0002114 150.8080  36.9326  1.00270880 14993"""
antenna_str = "Kuntunse, 5:45:2.48, -0:18:17.92, 116, 32.0"
raster_size_az = 15.0  # degrees
raster_size_el = 2.5  # degrees
cross_scan_extent_az = 6.0
cross_scan_extent_el = 1.0
number_of_scans = 75
scans_per_boresight = 5
scan_speed = 0.12
slew_speed = 0.2
settling_time = 15
boresight_time = 30
snake = True
plot = True
start_time_input = "2018-05-27 07:50:00"
output_filename = "scan6"


def slew(start_az, start_el, stop_az, stop_el, slew_speed, plot=False, labeltext="slew", plotformat="r."):
    """
    Slew from start_az/el to stop_az/el.

    Generate a linear set of az/el points from the start to the end, progressing according to a set speed.

    :param start_az: Starting azimuth, degrees.
    :type start_az: float
    :param start_el: Starting elevation, degrees.
    :type start_el: float
    :param stop_az: Ending azimuth, degrees.
    :type stop_az: float
    :param stop_el: Ending elevation, degrees.
    :type stop_el: float
    :param slew_speed: Absolute slew-speed (combination of azimuth and elevation according to Pythagoras),
                       degrees per second.
    :type slew_speed: float
    :param plot: Whether or not to plot the resulting pattern, default False.
    :type plot: bool
    :param labeltext: Custom text to add to the label for if a plot legend is used. Default is "slew".
    :type labeltext: str
    :param plotformat: Plot format to use for the slew if a plot is used. Default is "r."
    :type plotformat: str
    :return: Tuple containing azimuth points, elevation points, time for each point, and labels for each point.
    :rtype: np.ndarray(float), np.ndarray(float), np.ndarray(int), list

    """
    # Calculate distance, so that you can calculate the time.
    slew_distance = np.sqrt(np.square(stop_az - start_az) + np.square(stop_el - start_el))
    slew_time = int(np.floor(slew_distance / slew_speed)) + 1

    # Generate the az and el path from start to end. Nice and easy.
    az = np.linspace(start_az, stop_az, num=slew_time, endpoint=True, dtype=np.float)
    el = np.linspace(start_el, stop_el, num=slew_time, endpoint=True, dtype=np.float)

    # For the time being we're going to be issuing commands every second.
    # I'm of the opinion that it'll help make things smoother.
    t = np.ones(slew_time)

    # Making labeltext a parameter so that I can use the slew function for the scans.
    labels = [labeltext] * slew_time

    if plot:
        plt.plot(az, el, plotformat, label=labeltext)

    return az, el, t, labels


def dwell(dwell_az, dwell_el, dwell_time, plot=False, labeltext="dwell", plotformat="b*"):
    """
    Stay in one place for a while.

    Can be used for settling time or for boresight, with adjustable label text.

    :param dwell_az: Azimuth at which to dwell.
    :type dwell_az: float
    :param dwell_el: Elevation at which to dwell.
    :type dwell_el: float
    :param dwell_time: Time required for the dwell.
    :type dwell_time: int
    :param plot: Whether or not to plot the resulting pattern, default False.
    :type plot: bool
    :param labeltext: Custom text to add to the label for if a plot legend is used. Default is "slew".
    :type labeltext: str
    :param plotformat: Plot format to use for the slew if a plot is used. Default is "r."
    :type plotformat: str
    :return: Tuple containing azimuth points, elevation points, time for each point, and labels for each point.
    :rtype: np.ndarray(float), np.ndarray(float), np.ndarray(int), list
    """
    az = np.ones(dwell_time) * dwell_az
    el = np.ones(dwell_time) * dwell_el
    t = np.ones(dwell_time)
    labels = [labeltext] * dwell_time

    if plot:
        plt.plot(az, el, plotformat, label=labeltext)

    return az, el, t, labels


def generate_subscan(start_az, start_el, stop_az, stop_el, n_scans, scan_speed, slew_speed, settling_time, snake=True,
                     plot=False):
    """
    Generate subscan.

    A subscan is the number of back-and-forth scans between boresights.
    Assumption that the scan will always start at the bottom left and finish at the top, either left or right,
    depending on whether it "snakes" and the number of scans.

    :param start_az: Azimuth at bottom left corner of the block.
    :type start_az: float
    :param start_el: Elevation at the bottom left corner of the block.
    :type start_el: float
    :param stop_az: Azimuth at the top right corner of the block.
    :type stop_az: float
    :param stop_el: Elevation at the top right corner of the block.
    :type stop_el: float
    :param n_scans: Number of scans into which to divide the block.
    :type n_scans: int
    :param scan_speed: Rate at which to scan in azimuth, degrees per second.
    :type scan_speed: float
    :param slew_speed: Rate at which to slew in az / el (resultant), degrees per second.
    :type slew_speed: float
    :param settling_time: Amount of time to allocate for the dynamics of the telescope to settle, seconds.
    :type settling_time: int
    :param snake: Whether or not to go back and forth, i.e. not always the same direction to scan. Default True.
    :type snake: bool
    :param plot: Whether or not to plot the resulting pattern, default False.
    :type plot: bool
    :return: Tuple containing azimuth points, elevation points, time for each point, and labels for each point.
    :rtype: np.ndarray(float), np.ndarray(float), np.ndarray(int), list
    """
    az = np.array([])
    el = np.array([])
    t = np.array([])
    labels = []  # Can normal Python lists be concatenated too?

    cache = dwell(start_az, start_el,
                  settling_time, plot=plot, labeltext="settling")
    az = np.concatenate((az, cache[0]))
    el = np.concatenate((el, cache[1]))
    t = np.concatenate((t, cache[2]))
    labels = labels + cache[3]

    # Corner case.
    if n_scans == 1:
        cache = slew(start_az, start_el,
                     stop_az, start_el,
                     scan_speed, plot=plot, labeltext="scan", plotformat="g.")
        az = np.concatenate((az, cache[0]))
        el = np.concatenate((el, cache[1]))
        t = np.concatenate((t, cache[2]))
        labels = labels + cache[3]

    else:
        scan_spacing = (stop_el - start_el) / (n_scans - 1)

        for i in range(n_scans):
            if i % 2 == 0:
                cache = slew(start_az, start_el + i*scan_spacing,
                             stop_az, start_el + i*scan_spacing,
                             scan_speed, plot=plot, labeltext="scan", plotformat="g.")
                az = np.concatenate((az, cache[0]))
                el = np.concatenate((el, cache[1]))
                t = np.concatenate((t, cache[2]))
                labels = labels + cache[3]

            else:
                if snake:
                    cache = slew(stop_az, start_el + (i-1)*scan_spacing,
                                 stop_az, start_el + i*scan_spacing,
                                 slew_speed, plot=plot)  # This is a slew so default label stuff is fine.
                    az = np.concatenate((az, cache[0]))
                    el = np.concatenate((el, cache[1]))
                    t = np.concatenate((t, cache[2]))
                    labels = labels + cache[3]

                    cache = dwell(stop_az, start_el + i*scan_spacing,
                                  settling_time, plot=plot, labeltext="settling")
                    az = np.concatenate((az, cache[0]))
                    el = np.concatenate((el, cache[1]))
                    t = np.concatenate((t, cache[2]))
                    labels = labels + cache[3]

                    cache = slew(stop_az, start_el + i*scan_spacing,
                                 start_az, start_el + i*scan_spacing,
                                 scan_speed, plot=plot, labeltext="scan", plotformat="g.")
                    az = np.concatenate((az, cache[0]))
                    el = np.concatenate((el, cache[1]))
                    t = np.concatenate((t, cache[2]))
                    labels = labels + cache[3]

                    if i != (n_scans - 1):
                        cache = slew(start_az, start_el + i*scan_spacing,
                                     start_az, start_el + (i + 1)*scan_spacing,
                                     slew_speed, plot=plot)  # This is a slew so default label stuff is fine.
                        az = np.concatenate((az, cache[0]))
                        el = np.concatenate((el, cache[1]))
                        t = np.concatenate((t, cache[2]))
                        labels = labels + cache[3]

                        cache = dwell(start_az, start_el + (i + 1) * scan_spacing,
                                      settling_time, plot=plot, labeltext="settling")
                        az = np.concatenate((az, cache[0]))
                        el = np.concatenate((el, cache[1]))
                        t = np.concatenate((t, cache[2]))
                        labels = labels + cache[3]

                else:
                    cache = slew(stop_az, start_el + (i-1)*scan_spacing,
                                 start_az, start_el + i*scan_spacing,
                                 slew_speed, plot=plot)  # This is a slew so default label stuff is fine.
                    az = np.concatenate((az, cache[0]))
                    el = np.concatenate((el, cache[1]))
                    t = np.concatenate((t, cache[2]))
                    labels = labels + cache[3]

                    cache = dwell(start_az, start_el + i*scan_spacing,
                                  settling_time, plot=plot, labeltext="settling")
                    az = np.concatenate((az, cache[0]))
                    el = np.concatenate((el, cache[1]))
                    t = np.concatenate((t, cache[2]))
                    labels = labels + cache[3]

                    cache = slew(start_az, start_el + i*scan_spacing,
                                 stop_az, start_el + i*scan_spacing,
                                 scan_speed, plot=plot, labeltext="scan", plotformat="g.")
                    az = np.concatenate((az, cache[0]))
                    el = np.concatenate((el, cache[1]))
                    t = np.concatenate((t, cache[2]))
                    labels = labels + cache[3]

                    if i != (n_scans - 1):
                        cache = slew(stop_az, start_el + i*scan_spacing,
                                     start_az, start_el + (i + 1)*scan_spacing,
                                     slew_speed, plot=plot)  # This is a slew so default label stuff is fine.
                        az = np.concatenate((az, cache[0]))
                        el = np.concatenate((el, cache[1]))
                        t = np.concatenate((t, cache[2]))
                        labels = labels + cache[3]

                        cache = dwell(start_az, start_el + (i + 1)*scan_spacing,
                                      settling_time, plot=plot, labeltext="settling")
                        az = np.concatenate((az, cache[0]))
                        el = np.concatenate((el, cache[1]))
                        t = np.concatenate((t, cache[2]))
                        labels = labels + cache[3]

    cache = dwell(stop_az, stop_el,
                  settling_time, plot=plot, labeltext="settling")
    az = np.concatenate((az, cache[0]))
    el = np.concatenate((el, cache[1]))
    t = np.concatenate((t, cache[2]))
    labels = labels + cache[3]

    return az, el, t, labels


def generate_boresight(start_az, start_el, stop_az, stop_el, slew_speed, settling_time, boresight_time, plot=False):
    """
    Return to boresight before continuing.

    This function generates a slew to the origin (boresight), settles, dwells on the boresight for a specified time,
    then slews to the next point.

    :param start_az: Starting azimuth, degrees.
    :type start_az: float
    :param start_el: Starting elevation, degrees.
    :type start_el: float
    :param stop_az: Ending azimuth, degrees.
    :type stop_az: float
    :param stop_el: Ending elevation, degrees.
    :type stop_el: float
    :param slew_speed: Speed at which to slew, degrees per second.
    :type slew_speed: float
    :param settling_time: Time to allow the antenna to settle before starting the next thing, seconds.
    :type settling_time: int
    :param boresight_time: Time to stare at the boresight, seconds.
    :type boresight_time: int
    :param plot: Whether to add the scans to a plot.
    :type plot: bool
    :return: Tuple containing azimuth points, elevation points, time for each point, and labels for each point.
    :rtype: np.ndarray(float), np.ndarray(float), np.ndarray(int), list
    """
    az = np.array([])
    el = np.array([])
    t = np.array([])
    labels = []

    cache = slew(start_az, start_el, 0, 0, slew_speed, plot=plot)
    az = np.concatenate((az, cache[0]))
    el = np.concatenate((el, cache[1]))
    t = np.concatenate((t, cache[2]))
    labels = labels + cache[3]

    cache = dwell(0, 0, settling_time, plot=plot, labeltext="settling")
    az = np.concatenate((az, cache[0]))
    el = np.concatenate((el, cache[1]))
    t = np.concatenate((t, cache[2]))
    labels = labels + cache[3]

    cache = dwell(0, 0, boresight_time, plot=plot, labeltext="boresight")
    az = np.concatenate((az, cache[0]))
    el = np.concatenate((el, cache[1]))
    t = np.concatenate((t, cache[2]))
    labels = labels + cache[3]

    cache = slew(0, 0, stop_az, stop_el, slew_speed, plot=plot)
    az = np.concatenate((az, cache[0]))
    el = np.concatenate((el, cache[1]))
    t = np.concatenate((t, cache[2]))
    labels = labels + cache[3]

    return az, el, t, labels


def generate_holography_scan(az_extent, el_extent, num_scans, scans_per_boresight, scan_speed, slew_speed,
                             settling_time, boresight_time, snake=True, plot=False):
    """

    :param az_extent:
    :type az_extent:
    :param el_extent:
    :type el_extent:
    :param num_scans: Total number of scans into which to divide the extent.
    :type num_scans: int
    :param scans_per_boresight: Number of azimuth scans to do before returning to boresight.
    :type scans_per_boresight: int
    :param scan_speed: Speed at which to scan, degrees per second.
    :type scan_speed: float
    :param slew_speed: Speed at which to slew, degrees per second.
    :type slew_speed: float
    :param settling_time: Time to allow the antenna to settle before starting the next thing, seconds.
    :type settling_time: int
    :param boresight_time: Time to stare at the boresight, seconds.
    :type boresight_time: int
    :param snake:
    :type snake:
    :param plot: Whether to add the scans to a plot.
    :type plot: bool
    :return: Tuple containing azimuth points, elevation points, time for each point, and labels for each point.
    :rtype: np.ndarray(float), np.ndarray(float), np.ndarray(int), list
    """

    az = np.array([])
    el = np.array([])
    t = np.array([])
    labels = []

    scan_spacing = el_extent / (num_scans - 1)
    num_subscans = int(np.floor(num_scans / scans_per_boresight))
    subscan_sizing = scan_spacing * (scans_per_boresight - 1)
    partial_subscan_size = num_scans % scans_per_boresight

    # This is so that the boresight slew comes from the right place at the end of a subscan.
    sign = 1
    if scans_per_boresight % 2 == 0:
        sign = -1

    # Corner case:
    if scans_per_boresight == 1:
        raise NotImplementedError("Haven't gotten to this corner case yet.")
    # TODO: implement this. Shouldn't be too hard...

    for i in range(num_subscans):
        cache = generate_subscan(-az_extent/2, -el_extent/2 + i*(subscan_sizing + scan_spacing),
                         az_extent/2,  -el_extent/2 + (i + 1)*(subscan_sizing + scan_spacing) - scan_spacing,
                         scans_per_boresight, scan_speed, slew_speed, settling_time, snake=snake, plot=plot)
        az = np.concatenate((az, cache[0]))
        el = np.concatenate((el, cache[1]))
        t = np.concatenate((t, cache[2]))
        labels = labels + cache[3]

        if not (i == num_subscans - 1 and partial_subscan_size == 0):
            cache = generate_boresight(sign * az_extent/2,  -el_extent/2 + (i + 1)*(subscan_sizing + scan_spacing) - scan_spacing,
                                       -az_extent / 2, -el_extent / 2 + (i + 1) * (subscan_sizing + scan_spacing),
                                       slew_speed, settling_time, boresight_time, plot=plot)
            az = np.concatenate((az, cache[0]))
            el = np.concatenate((el, cache[1]))
            t = np.concatenate((t, cache[2]))
            labels = labels + cache[3]

    if partial_subscan_size != 0:
        cache = generate_subscan(-az_extent/2, el_extent/2 - (partial_subscan_size - 1)*(scan_spacing),
                                 az_extent/2, el_extent/2, partial_subscan_size, scan_speed, slew_speed, settling_time,
                                 snake=snake, plot=plot)
        az = np.concatenate((az, cache[0]))
        el = np.concatenate((el, cache[1]))
        t = np.concatenate((t, cache[2]))
        labels = labels + cache[3]

    return az, el, t, labels


def generate_cross_scan(start_az, start_el, stop_az, stop_el, scan_extent_az, scan_extent_el, scan_speed, slew_speed,
                        settling_time, boresight_time, plot=False):
    """
    Scan across the target.

    Detects whether you're starting at the origin, and if not, slews there first, settles, then starts the cross-scan.

    :param start_az: Starting azimuth, degrees.
    :type start_az: float
    :param start_el: Starting elevation, degrees.
    :type start_el: float
    :param stop_az: Ending azimuth, degrees.
    :type stop_az: float
    :param stop_el: Ending elevation, degrees.
    :type stop_el: float
    :param scan_extent_az: The size of the scan in azimuth, degrees.
    :type scan_extent_az: float
    :param scan_extent_el: The size of the scan in elevation, degrees.
    :type scan_extent_el: float
    :param scan_speed: Speed at which to scan, degrees per second.
    :type scan_speed: float
    :param slew_speed: Speed at which to slew, degrees per second.
    :type slew_speed: float
    :param settling_time: Time to allow the antenna to settle before starting the next thing, seconds.
    :type settling_time: int
    :param boresight_time: Time to stare at the boresight, seconds.
    :type boresight_time: int
    :param plot: Whether to add the scans to a plot.
    :type plot: bool
    :return: Tuple containing azimuth points, elevation points, time for each point, and labels for each point.
    :rtype: np.ndarray(float), np.ndarray(float), np.ndarray(int), list
    """
    az = np.array([])
    el = np.array([])
    t = np.array([])
    labels = []

    # If we start away from boresight, go there.
    if (start_az, start_el) != (0, 0):
        cache = slew(start_az, start_el,
                     0, 0,
                     slew_speed, plot=plot)
        az = np.concatenate((az, cache[0]))
        el = np.concatenate((el, cache[1]))
        t = np.concatenate((t, cache[2]))
        labels = labels + cache[3]

        cache = dwell(0, 0, settling_time, plot=plot, labeltext="settling")
        az = np.concatenate((az, cache[0]))
        el = np.concatenate((el, cache[1]))
        t = np.concatenate((t, cache[2]))
        labels = labels + cache[3]

    # Stare at boresight for a bit.
    cache = dwell(0, 0, boresight_time, plot=plot, labeltext="boresight")
    az = np.concatenate((az, cache[0]))
    el = np.concatenate((el, cache[1]))
    t = np.concatenate((t, cache[2]))
    labels = labels + cache[3]

    # Slew to the start point, top of elev extent.
    cache = slew(0, 0,
                 0, scan_extent_el/2,
                 slew_speed, plot=plot)
    az = np.concatenate((az, cache[0]))
    el = np.concatenate((el, cache[1]))
    t = np.concatenate((t, cache[2]))
    labels = labels + cache[3]

    # Settle.
    cache = dwell(0, scan_extent_el/2, settling_time, plot=plot, labeltext="settling")
    az = np.concatenate((az, cache[0]))
    el = np.concatenate((el, cache[1]))
    t = np.concatenate((t, cache[2]))
    labels = labels + cache[3]

    # Scan through to the bottom.
    cache = slew(0, scan_extent_el/2,
                 0, -scan_extent_el / 2,
                 scan_speed, plot=plot, labeltext="cross-scan", plotformat="g.")
    az = np.concatenate((az, cache[0]))
    el = np.concatenate((el, cache[1]))
    t = np.concatenate((t, cache[2]))
    labels = labels + cache[3]

    # Settle.
    cache = dwell(0, -scan_extent_el / 2, settling_time, plot, labeltext="settling")
    az = np.concatenate((az, cache[0]))
    el = np.concatenate((el, cache[1]))
    t = np.concatenate((t, cache[2]))
    labels = labels + cache[3]

    # Go to the azimuth starting point.
    cache = slew(0, -scan_extent_el / 2,
                 -scan_extent_az/2, 0,
                 slew_speed, plot)
    az = np.concatenate((az, cache[0]))
    el = np.concatenate((el, cache[1]))
    t = np.concatenate((t, cache[2]))
    labels = labels + cache[3]

    # Settle.
    cache = dwell(-scan_extent_az/2, 0, settling_time, plot, labeltext="settlling")
    az = np.concatenate((az, cache[0]))
    el = np.concatenate((el, cache[1]))
    t = np.concatenate((t, cache[2]))
    labels = labels + cache[3]

    # Scan through to the other side.
    cache = slew(-scan_extent_az/2, 0,
                 scan_extent_az/2, 0,
                 scan_speed, plot, labeltext="cross-scan", plotformat="g.")
    az = np.concatenate((az, cache[0]))
    el = np.concatenate((el, cache[1]))
    t = np.concatenate((t, cache[2]))
    labels = labels + cache[3]

    # Settle.
    cache = dwell(scan_extent_az/2, 0, settling_time, plot, labeltext="settling")
    az = np.concatenate((az, cache[0]))
    el = np.concatenate((el, cache[1]))
    t = np.concatenate((t, cache[2]))
    labels = labels + cache[3]

    # Go to the ending point.
    cache = slew(scan_extent_az/2, 0,
                 stop_az, stop_el,
                 slew_speed, plot)
    az = np.concatenate((az, cache[0]))
    el = np.concatenate((el, cache[1]))
    t = np.concatenate((t, cache[2]))
    labels = labels + cache[3]

    return az, el, t, labels


if __name__ == "__main__":

    start_time = 0
    try:
        start_time = calendar.timegm(time.strptime(start_time_input, "%Y-%m-%d %H:%M:%S"))
    except ValueError:
        print "Either start date or time are incorrectly formatted. Use format: 2001-09-11 09:57:55"
        exit(-1)

    myTarget = katpoint.Target("%s, tle, %s" % (satellite_name, TLE))
    antenna = katpoint.Antenna(antenna_str)
    myTarget.antenna = antenna

    az = np.array([])
    el = np.array([])
    t = np.array([])
    labels = []

    cache = generate_cross_scan(0, 0,
                                -raster_size_az/2, -raster_size_el/2,
                                cross_scan_extent_az, cross_scan_extent_el,
                                scan_speed, slew_speed, settling_time, boresight_time, plot=plot)
    az = np.concatenate((az, cache[0]))
    el = np.concatenate((el, cache[1]))
    t = np.concatenate((t, cache[2]))
    labels = labels + cache[3]

    cache = generate_holography_scan(raster_size_az, raster_size_el,
                                     number_of_scans, scans_per_boresight,
                                     scan_speed, slew_speed, settling_time,
                                     boresight_time, snake, plot)
    az = np.concatenate((az, cache[0]))
    el = np.concatenate((el, cache[1]))
    t = np.concatenate((t, cache[2]))
    labels = labels + cache[3]

    cache = generate_cross_scan(raster_size_az/2, raster_size_el/2,
                                0, 0,
                                cross_scan_extent_az, cross_scan_extent_el,
                                scan_speed, slew_speed, settling_time, boresight_time, plot=plot)
    az = np.concatenate((az, cache[0]))
    el = np.concatenate((el, cache[1]))
    t = np.concatenate((t, cache[2]))
    labels = labels + cache[3]

    # Calculate the time that we want to be at each given point.
    t = np.cumsum(t)
    t -= t[0]
    t += start_time

    # Find where in real-space the satellite is.
    target_position_az = []
    target_position_el = []
    for time_step in t:
        target_position_az.append(np.degrees(myTarget.azel(time_step)[0]))
        target_position_el.append(np.degrees(myTarget.azel(time_step)[1]))
    target_position_az = np.array(target_position_az)
    target_position_el = np.array(target_position_el)

    # Add the relative position values to the absolute position of the satellite.
    az += target_position_az
    el += target_position_el

    if plot:
        plt.grid()
        plt.show()

        plt.figure()
        plt.plot(target_position_az, target_position_el)
        plt.show()

    with open("%s.snp" % output_filename, "w") as output_file:
        output_file.write("\" Holography raster script developed 2018-05-25 by James Smith.\n")
        output_file.write("\" Field system snap file automagically generated at: %s.\n"
                          % time.strftime("%Y-%m-%d %H:%M:%S"))
        output_file.write("\" File: %s.snp\n" % output_filename)
        output_file.write("\" Target: %s\n" % repr(myTarget))
        output_file.write("\" Target nominal az/el: %.2f, %.2f\n" % (np.average(target_position_az),
                                                                     np.average(target_position_el)))
        output_file.write("\" Raster size: %.2f x %.2f degrees\n" % (raster_size_az, raster_size_el))
        output_file.write("\" Schedule to start at: %s\n" % start_time_input)
        # TODO: Consider writing a "finish time" field as well.
        for i in range(len(t)):
            output_file.write("\"#%d,%s,nominal\n!%s\nsource=azel,%.6fd,%.6fd\n" %
                              (t[i], labels[i],
                               time.strftime("%Y.%j.%H:%M:%S", time.gmtime(t[i])),
                               az[i], el[i]))
        output_file.write("stop")
