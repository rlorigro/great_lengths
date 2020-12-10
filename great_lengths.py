from modules.plot import plot_iterative_histogram, plot_ngx
from modules.IterativeHistogram import IterativeHistogram
from collections import Counter
from subprocess import run
import argparse
import sys
import os


def write_report(quartiles, output_dir):
    path = os.path.join(output_dir, "report.tsv")
    sys.stderr.write("SAVING REPORT: %s\n" % path)

    with open(path, 'w') as output_file:
        line = "quartiles" + '\t' + '\t'.join(list(map(str, quartiles)))
        output_file.write(line)


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

            if i_quartile == len(quartiles):
                break

        n_visited += frequency

    print(quartile_values)

    # Append the min and max
    quartile_values = [length_frequencies[0][0]] + quartile_values + [length_frequencies[-1][0]]

    return quartile_values


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
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

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

    if histogram_max == histogram_min:
        exit("ERROR: cannot create histogram for read distribution containing only one length")

    histogram = IterativeHistogram(start=histogram_min, stop=histogram_max, n_bins=histogram_n_bins)

    length_frequencies = Counter()
    total_length = 0
    n_items = 0
    with open(index_path) as file:
        for line in file:
            print(length)
            length = int(line.split('\t')[1])
            length_frequencies[length] += 1
            histogram.update(length)
            total_length += length
            n_items += 1

    length_frequencies = sorted(length_frequencies.items(), key=lambda x: x[0])

    quartiles = find_quartiles(length_frequencies, n_items=n_items)

    plot_iterative_histogram(histogram, output_dir=output_dir)
    plot_ngx(length_frequencies, total_length=total_length, output_dir=output_dir)
    write_report(quartiles=quartiles, output_dir=output_dir)



def string_as_bool(s):
    s = s.lower()
    boolean = None

    if s in {"t", "true", "1", "y", "yes"}:
        boolean = True
    elif s in {"f", "false", "0", "n", "no"}:
        boolean = False
    else:
        exit("Error: invalid argument specified for boolean flag: %s" % s)

    return boolean


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.register("type", "bool", string_as_bool)  # add type keyword to registries

    parser.add_argument(
        "--input","-i",
        type=str,
        required=True,
        help="path of file containing FASTA/FASTQ sequence file"
    )
    parser.add_argument(
        "--output_dir","-o",
        type=str,
        required=True,
        help="path of destination directory for output files"
    )
    parser.add_argument(
        "--hist_min",
        type=int,
        required=False,
        default=0,
        help="Minimum of histogram range (default=100,000)"
    )
    parser.add_argument(
        "--hist_max",
        type=int,
        required=False,
        default=100_000,
        help="Maximum of histogram range, (default=0)"
    )
    parser.add_argument(
        "--hist_n_bins",
        type=int,
        required=False,
        default=500,
        help="Maximum of histogram range, (default=1000)"
    )
    parser.add_argument(
        "--hist_auto_bounds",
        type=int,
        required=False,
        default=500,
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
