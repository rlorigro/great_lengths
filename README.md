# great_lengths
Very simple tool for generating read length distributions from fastq/fasta files

## Output
A histogram of read lengths (with automatic or predefined bounds):

![histogram](/images/histogram.svg)


A Nx plot of read lengths:

![Nx](/images/Nx.svg)


A summary table with statistics:
```
file	GM24385_1_Guppy_4.2.2_prom.fastq
total_reads	2900987
total_bp	26945078155
total_Gbp	26
min	1
max	503494
mean	9288
quartile_25	17
quartile_50	932
quartile_75	6555
N25	21001
N50	52031
N75	88313
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
