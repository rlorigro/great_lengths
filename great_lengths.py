from modules.plot import plot_iterative_histogram, plot_ngx
from modules.IterativeHistogram import IterativeHistogram
from collections import Counter
from subprocess import run
import argparse
import sys
import os



def find_quartiles(length_frequencies, n_items):
    n_visited = 0
    i_quartile = 0
    quartiles = [0.25, 0.5, 0.75]
    quartile_values = list()

    for length,frequency in length_frequencies:
        # This instance of 'length' may have occured multiple times, as measured by 'frequency'
        # Check if all the instances of this length would exceed any quartile
        while float(n_visited + frequency)/float(n_items) > quartiles[i_quartile]:
            quartile_values.append(length)
            i_quartile += 1

            # terminate loop
            if i_quartile == len(quartiles):
                break

        n_visited += frequency

        # exit
        if len(quartile_values) == len(quartiles):
            break

    # Append the min and max
    quartile_values = [length_frequencies[0][0]] + quartile_values + [length_frequencies[-1][0]]

    return quartile_values


def find_n25_n50_n75(sorted_length_frequencies, total_size):
    n25_threshold = 0.25 * total_size
    n50_threshold = 0.5 * total_size
    n75_threshold = 0.75 * total_size
    ns = list()
    current_total = 0

    assert(sorted_length_frequencies[0][0] <= sorted_length_frequencies[-1][0])
    for length,frequency in sorted_length_frequencies:
        current_total += length * frequency
        if len(ns) == 0 and current_total >= n25_threshold:
            ns.append(length)
        if len(ns) == 1 and current_total >= n50_threshold:
            ns.append(length)
        if len(ns) == 2 and current_total >= n75_threshold:
            ns.append(length)
            break
        assert(len(ns) < 3)
    assert(len(ns) == 3)

    return ns


# Use a system call to samtools faidx to build the index
def build_index(path):
    index_path = path + ".fai"
    if os.path.exists(index_path):
        sys.stderr.write("Index exists\n")
        if os.path.getmtime(index_path) < os.path.getmtime(path):
            sys.stderr.write("WARNING: fasta/q index path is older than file itself, may be out of date: " + index_path + "\n")
    else:
        sys.stderr.write("No index found, indexing... ")
        arguments = ["samtools", "faidx", path]
        run(arguments, check=True)
        sys.stderr.write("Done\n")

    sys.stderr.flush()

    return index_path


def main(input_path, output_dir, histogram_min, histogram_max, histogram_n_bins, use_auto_bounds):

    # sanity check
    if histogram_max == histogram_min:
        exit("ERROR: cannot create histogram for read distribution containing only one length")

    # ensure output dir if set
    if output_dir is not None:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
    else:
        sys.stderr.write("No output directory specified, will print report to stdout and not produce plots")

    # get index
    index_path = build_index(input_path)
    if use_auto_bounds:
        sys.stderr.write("WARNING: using auto bounds increases run time\n")
        histogram_min = sys.maxsize
        histogram_max = 0
        with open(index_path) as file:
            for line in file:
                length = int(line.split('\t')[1])

                if length > histogram_max:
                    histogram_max = length
                if length < histogram_min:
                    histogram_min = length
        sys.stderr.write("Using automatically determined bounds for histogram: [%d,%d]\n"%(histogram_min, histogram_max))


    histogram = IterativeHistogram(start=histogram_min, stop=histogram_max, n_bins=histogram_n_bins)

    length_frequencies = Counter()
    total_length = 0
    n_items = 0
    with open(index_path) as file:
        for line in file:
            length = int(line.split('\t')[1])
            length_frequencies[length] += 1
            histogram.update(length)
            total_length += length
            n_items += 1

    # get data
    length_frequencies = sorted(length_frequencies.items(), key=lambda x: x[0])
    quartiles = find_quartiles(length_frequencies, n_items=n_items)
    ns = find_n25_n50_n75(length_frequencies, total_length)

    # plots
    if output_dir is not None:
        plot_iterative_histogram(histogram, output_dir=output_dir)
        plot_ngx(length_frequencies, total_length=total_length, output_dir=output_dir)

    # write report
    output_file = None
    try:
        if output_dir is not None:
            path = os.path.join(output_dir, "report.tsv")
            sys.stderr.write("SAVING REPORT: %s\n" % path)
            output_file = open(path, 'w')
        else:
            output_file = sys.stdout

        print("file\t{}".format(os.path.basename(input_path)), file=output_file)
        print("total_reads\t{}".format(n_items), file=output_file)
        print("total_bp\t{}".format(total_length), file=output_file)
        print("total_Gbp\t{}".format(total_length // 1000000000), file=output_file)
        print("min\t{}".format(quartiles[0]), file=output_file)
        print("max\t{}".format(quartiles[4]), file=output_file)
        print("mean\t{}".format(total_length // n_items), file=output_file)
        print("quartile_25\t{}".format(quartiles[1]), file=output_file)
        print("quartile_50\t{}".format(quartiles[2]), file=output_file)
        print("quartile_75\t{}".format(quartiles[3]), file=output_file)
        print("N25\t{}".format(ns[0]), file=output_file)
        print("N50\t{}".format(ns[1]), file=output_file)
        print("N75\t{}".format(ns[2]), file=output_file)

    finally:
        if output_dir is not None and output_file is not None:
            output_file.close()



if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input","-i",
        type=str,
        required=True,
        help="path of file containing FASTA/FASTQ sequence file"
    )
    parser.add_argument(
        "--output_dir","-o",
        type=str,
        required=False,
        default=None,
        help="path of destination directory for output files"
    )
    parser.add_argument(
        "--hist_min",
        type=int,
        required=False,
        default=0,
        help="Minimum of histogram range (default=0)"
    )
    parser.add_argument(
        "--hist_max",
        type=int,
        required=False,
        default=100_000,
        help="Maximum of histogram range, (default=100000)"
    )
    parser.add_argument(
        "--hist_n_bins",
        type=int,
        required=False,
        default=500,
        help="Maximum of histogram range, (default=500)"
    )
    parser.add_argument(
        "--hist_auto_bounds",
        required=False,
        default=False,
        action='store_true',
        help="Automatically determine the histogram min/max bounds from the min/max in the data"
    )

    args = parser.parse_args()

    main(
        input_path=args.input,
        output_dir=args.output_dir,
        histogram_min=args.hist_min,
        histogram_max=args.hist_max,
        histogram_n_bins=args.hist_n_bins,
        use_auto_bounds=args.hist_auto_bounds
    )
