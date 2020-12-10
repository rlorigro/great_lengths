from modules.IterativeHistogram import IterativeHistogram
from modules.plot import plot_iterative_histogram
from pyfaidx import Faidx
import argparse
import sys
import os


def build_index(path):
    indexer = Faidx(path)

    index_path = path + ".fai"
    if os.path.exists(index_path):
        sys.stderr.write("Index exists\n")
        if os.path.getmtime(index_path) < os.path.getmtime(path):
            sys.stderr.write("WARNING: fasta/q index path is older than file itself, may be out of date: " + index_path + "\n")
    else:
        sys.stderr.write("No index found, indexing... ")
        indexer = indexer.build_index()
        sys.stderr.write("Done\n")

    sys.stderr.flush()

    return indexer


def main(input_path, output_dir, histogram_min, histogram_max, histogram_n_bins, use_auto_bounds):
    filename_prefix = (os.path.basename(input_path)).split('.')[0]

    indexer = build_index(input_path)
    lengths = list()

    if use_auto_bounds:
        # For some godforsaken reason this library has no min/max length method for FASTA
        sys.stderr.write("WARNING: using auto bounds can increase run time\n")

        histogram_min = sys.maxsize
        histogram_max = 0

        for index in indexer:
            if len(s) > histogram_max:
                histogram_max = len(s)
            if len(s) < histogram_min:
                histogram_min = len(s)

        sys.stderr.write("Using automatically determined bounds for histogram: [%d,%d]\n"%(histogram_min, histogram_max))

    if histogram_max == histogram_min:
        exit("ERROR: cannot create histogram for read distribution containing only one length")

    histogram = IterativeHistogram(start=histogram_min, stop=histogram_max, n_bins=histogram_n_bins)

    for sequence in file_parser:
        print(sequence.name, len(sequence))
        lengths.append(len(sequence))
        histogram.update(len(sequence))

    histogram_path_name = os.path.join(output_dir, filename_prefix + ".png")
    plot_iterative_histogram(histogram, path=histogram_path_name)


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
