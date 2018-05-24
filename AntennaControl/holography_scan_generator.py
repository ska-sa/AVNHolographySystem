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
raster_size_az = 1.5  # degrees
raster_size_el = 1.5  # degrees
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
        plt.plot(az, el, plotformat, label="slew")

    return az, el, t, labels


def generate_subscan(start_az, start_el, stop_az, stop_el, n_scans, scan_speed, slew_speed, settling_time, snake=True):
    # This is a place-holder.
    az = np.array([start_az, stop_az])
    el = np.array([start_el, stop_el])
    t = [1, 2]
    labels = ["foo", "bar"]
    plt.plot(az, el)

    return az, el, t, labels


if __name__ == "__main__":
    plt.figure()
    generate_subscan(0, 0, 10, 2, 3, 1, 1, 1)
    bar = slew(10, 2, 0, 3, 1, plot=True)
    print bar
    generate_subscan(0, 3, 10, 5, 3, 1, 1, 1)

    plt.show()
