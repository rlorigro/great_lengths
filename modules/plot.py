from matplotlib import pyplot
import numpy
import sys
import os


def plot_ngx(lengths, output_dir):
    total_length = sum(lengths)

    figure = pyplot.figure()
    axes = pyplot.axes()

    legend_names = list()

    x1 = 0
    y_prev = None

    x_coords = list()
    y_coords = list()

    for length in lengths:
        y = length
        width = float(length) / float(total_length)
        x2 = x1 + width

        if y_prev is not None:
            x_coords.extend([x1, x1])
            y_coords.extend([y_prev, y])

        x_coords.extend([x1, x2])
        y_coords.extend([y, y])

        x1 = x2
        y_prev = y

    if y_coords[-1] != 0:
        y_coords.append(0)
        x_coords.append(x_coords[-1])

    axes.plot(x_coords, y_coords, linewidth=0.6)

    axes.legend(legend_names)

    axes.axvline(0.5, linestyle="--", alpha=0.3, linewidth=0.7, zorder=-1)

    axes.set_xlim([0, 1])

    axes.set_title("NGx")
    axes.set_xlabel("Cumulative Coverage")

    print("SAVING FIGURE: %s" % path)
    figure.savefig(path + ".png", dpi=300)
    figure.savefig(path + ".pdf", dpi=300)

    pyplot.close()

def plot_iterative_histogram(iterative_histogram, path=None):
    figure = pyplot.figure()
    axes = pyplot.axes()

    bounds = numpy.array(iterative_histogram.edges)

    center = (bounds[:-1] + bounds[1:]) / 2

    axes.bar(center, iterative_histogram.histogram, width=iterative_histogram.bin_size, align="center")

    axes.set_xlabel("Read length (bp)")
    axes.set_ylabel("Frequency")

    if path is not None:
        figure.savefig(path, dpi=200)

    pyplot.close()
