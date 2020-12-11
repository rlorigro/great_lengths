# great_lengths
Very simple tool for generating read length distributions from fastq/fasta files

## Output
A histogram of read lengths (with automatic or predefined bounds):

![histogram](/images/histogram.svg)


A Nx plot of read lengths:

![Nx](/images/Nx.svg)


A summary table with statistics:
```
file    simple.fastq
total_reads     22
total_bp        1042
total_Gbp       0
min     12
max     76
mean    47
quartile_25     12
quartile_50     76
quartile_75     76
N25     76
N50     76
N75     76
```

An unabridged table of all observed lengths:
```
Length  Frequency
12      9
22      1
76      12
```

## Usage
```
usage: great_lengths.py [-h] --input INPUT [--output_dir OUTPUT_DIR]
                        [--hist_min HIST_MIN] [--hist_max HIST_MAX]
                        [--hist_n_bins HIST_N_BINS] [--hist_auto_bounds]
                        [--unabridged]

optional arguments:
  -h, --help            show this help message and exit
  --input INPUT, -i INPUT
                        path of file containing FASTA/FASTQ sequence file
  --output_dir OUTPUT_DIR, -o OUTPUT_DIR
                        path of destination directory for output files
  --hist_min HIST_MIN   Minimum of histogram range (default=0)
  --hist_max HIST_MAX   Maximum of histogram range, (default=100000)
  --hist_n_bins HIST_N_BINS
                        Maximum of histogram range, (default=500)
  --hist_auto_bounds    Automatically determine the histogram min/max bounds
                        from the min/max in the data
  --unabridged, -u      Dump the full distribution of all observed lengths in
                        tsv with columns: length,frequency
```

## Dependencies
The main script expects samtools 1.9 or greater to be installed and added to the PATH environment variable. Version 1.7 will fail to index FASTQ files and should not be used.

Python 3
 - numpy
 - matplotlib
