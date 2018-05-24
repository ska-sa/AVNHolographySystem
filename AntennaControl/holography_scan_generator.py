import numpy as np
import matplotlib.pyplot as plt
import time
import calendar  # for the one timegm function...
import katpoint

# Edit these before each run.
satellite_name = "SES-5"
TLE = """SES-5
1 38652C 12036A   18142.65625000  .00000124  00000-0  00000-0 0  1423
2 38652   0.0271 263.9528 0001999 149.8080  67.7044  1.00269513    15"""
antenna_str = "Kuntunse, 5:45:2.48, -0:18:17.92, 116, 32.0"
raster_size_az = 10.0  # degrees
raster_size_el = 10.0  # degrees
number_of_scans = 30
scans_per_boresight = 5
scan_speed = 0.01
slew_speed = 0.01
settling_time = 5
boresight_time = 30
snake = True
plot = True
start_time_input = "2018-05-23 20:30:00"
output_filename = "output"


def slew(start_az, start_el, stop_az, stop_el, slew_speed, plot=False, labeltext="slew", plotformat="r."):
    """
    Slew from start_az/el to stop_az/el.

    Generate a linear set of az/el points from the start to the end, progressing according to a set speed.

    :param start_az: Azimuth coordinate of the starting point, degrees clockwise of N.
    :type start_az: float
    :param start_el: Elevation coordinate of the starting point, degrees upwards from horizon.
    :type start_el: float
    :param stop_az: Azimuth coordinate of the finishing point, degrees clockwise of N.
    :type stop_az: float
    :param stop_el: Elevation coordinate of the finishing point, degrees upwards from horizon.
    :type stop_el: float
    :param slew_speed: Absolute slew-speed (combination of azimuth and elevation according to Pythagoras),
                       degrees per second.
    :type slew_speed: float
    :param plot: Add the slew to a plot if True.
    :type plot: bool
    :param labeltext: Custom text to add to the label for if a plot legend is used. Default is "slew".
    :type labeltext: str
    :param plotformat: Plot format to use for the slew if a plot is used. Default is "r."
    :type plotformat: str
    :return az
    :return el

    """
    # TODO finish the docstring's return info.
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

    :param dwell_az:
    :type dwell_az:
    :param dwell_el:
    :type dwell_el:
    :param dwell_time:
    :type dwell_time:
    :param plot:
    :type plot:
    :param labeltext:
    :type labeltext:
    :param plotformat:
    :type plotformat:
    :return:
    :rtype:
    """
    # TODO: Finish docstring.
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

    # Corner case.
    if n_scans == 1:
        cache = slew(start_az, start_el,
                     stop_az, start_el,
                     scan_speed, plot=plot, labeltext="scan", plotformat="g.")
        az = np.concatenate((az, cache[0]))
        el = np.concatenate((el, cache[1]))
        t = np.concatenate((t, cache[2]))
        labels = labels + cache[3]

    # TODO: Continue to handle these return values.

    else:
        scan_spacing = (stop_el - start_el) / (n_scans - 1)

        for i in range(n_scans):
            if i % 2 == 0:
                slew(start_az, start_el + i*scan_spacing,
                     stop_az, start_el + i*scan_spacing,
                     scan_speed, plot=plot, labeltext="scan", plotformat="g.")
            else:
                if snake:
                    slew(stop_az, start_el + (i-1)*scan_spacing,
                         stop_az, start_el + i*scan_spacing,
                         slew_speed, plot=plot)  # This is a slew so default label stuff is fine.
                    dwell(stop_az, start_el + i*scan_spacing,
                          settling_time, plot=plot, labeltext="settling")
                    slew(stop_az, start_el + i*scan_spacing,
                         start_az, start_el + i*scan_spacing,
                         scan_speed, plot=plot, labeltext="scan", plotformat="g.")
                    if i != (n_scans - 1):
                        slew(start_az, start_el + i*scan_spacing,
                             start_az, start_el + (i + 1)*scan_spacing,
                             slew_speed, plot=plot)  # This is a slew so default label stuff is fine.
                        dwell(start_az, start_el + (i + 1) * scan_spacing,
                              settling_time, plot=plot, labeltext="settling")
                else:
                    slew(stop_az, start_el + (i-1)*scan_spacing,
                         start_az, start_el + i*scan_spacing,
                         slew_speed, plot=plot)  # This is a slew so default label stuff is fine.
                    dwell(start_az, start_el + i*scan_spacing,
                          settling_time, plot=plot, labeltext="settling")
                    slew(start_az, start_el + i*scan_spacing,
                         stop_az, start_el + i*scan_spacing,
                         scan_speed, plot=plot, labeltext="scan", plotformat="g.")
                    if i != (n_scans - 1):
                        slew(stop_az, start_el + i*scan_spacing,
                             start_az, start_el + (i + 1)*scan_spacing,
                             slew_speed, plot=plot)  # This is a slew so default label stuff is fine.
                        dwell(start_az, start_el + (i + 1)*scan_spacing,
                              settling_time, plot=plot, labeltext="settling")

    return az, el, t, labels


def generate_boresight(start_az, start_el, stop_az, stop_el, slew_speed, settling_time, dwell_time, plot=False):
    """

    :param start_az:
    :type start_az:
    :param start_el:
    :type start_el:
    :param stop_az:
    :type stop_az:
    :param stop_el:
    :type stop_el:
    :param slew_speed:
    :type slew_speed:
    :param settling_time:
    :type settling_time:
    :param dwell_time:
    :type dwell_time:
    :param plot:
    :type plot:
    :return:
    :rtype:
    """
    # TODO: finish docstring

    az = np.array([])
    el = np.array([])
    t = np.array([])
    labels = []  # Can normal Python lists be concatenated too?

    slew(start_az, start_el, 0, 0, slew_speed, plot=plot)
    dwell(0, 0, settling_time, plot=plot, labeltext="settling")
    dwell(0, 0, dwell_time, plot=plot, labeltext="boresight")
    slew(0, 0, stop_az, stop_el, slew_speed, plot=plot)

    return az, el, t, labels


def generate_holography_scan(az_extent, el_extent, num_scans, scans_per_boresight, scan_speed, slew_speed,
                             settling_time, boresight_time, snake=True, plot=False):
    az = np.array([])
    el = np.array([])
    t = np.array([])
    labels = []  # Can normal Python lists be concatenated too?

    scan_spacing = el_extent / (num_scans - 1)
    num_subscans = int(np.floor(num_scans / scans_per_boresight))
    subscan_sizing = scan_spacing * (scans_per_boresight - 1)
    partial_subscan_size = num_scans % scans_per_boresight

    sign = 1
    if scans_per_boresight % 2 == 0:
        sign = -1

    # TODO: Insert a cross scan here.

    # Corner case:
    if scans_per_boresight == 1:
        raise NotImplementedError("Haven't gotten to this corner case yet.")
      # TODO: implement this. Shouldn't be too hard...


    i = 0
    for i in range(num_subscans):
        generate_subscan(-az_extent/2, -el_extent/2 + i*(subscan_sizing + scan_spacing),
                         az_extent/2,  -el_extent/2 + (i + 1)*(subscan_sizing + scan_spacing) - scan_spacing,
                         scans_per_boresight, scan_speed, slew_speed, settling_time, snake=snake, plot=plot)

        if not (i == num_subscans - 1 and partial_subscan_size == 0):
            generate_boresight(sign * az_extent/2,  -el_extent/2 + (i + 1)*(subscan_sizing + scan_spacing) - scan_spacing,
                               -az_extent / 2, -el_extent / 2 + (i + 1) * (subscan_sizing + scan_spacing),
                               slew_speed, settling_time, boresight_time, plot=plot)

    if partial_subscan_size != 0:
        generate_subscan(-az_extent/2, el_extent/2 - (partial_subscan_size - 1)*(scan_spacing),
                         az_extent/2, el_extent/2, partial_subscan_size, scan_speed, slew_speed, settling_time,
                         snake=snake, plot=plot)


    # TODO: Insert another cross scan here.

    return az, el, t, labels


if __name__ == "__main__":

    d_az, d_el, d_time, labels = generate_holography_scan(raster_size_az, raster_size_el,
                                                          number_of_scans, scans_per_boresight,
                                                          scan_speed, slew_speed, settling_time,
                                                          boresight_time, snake, plot)

    # TODO: Write out the file, as previously. Shouldn't change much.

    if plot:
        plt.grid()
        plt.show()
