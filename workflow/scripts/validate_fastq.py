#!/usr/bin/env python3

import gzip
import re
import sys
from pathlib import Path


def count_gzip_lines(filepath):
    """Count the number of lines in a gzipped file."""
    with gzip.open(filepath, 'rt') as f:
        return sum(1 for _ in f)


def normalize_read_name(read_name):
    """Normalize read name by removing '/1' or '/2' suffixes."""
    return re.sub(r'(/[12])$', '', read_name)


def validate_fastq_single(input_file, output_file):
    """Validate single-end FASTQ file."""
    n0 = count_gzip_lines(input_file)
    if n0 % 4 == 0:
        Path(output_file).touch()
    else:
        print(f"ERROR: Number of lines in FASTQ file is not a multiple of four: {n0}", file=sys.stderr)
        sys.exit(1)


def validate_fastq_paired(input_file1, input_file2, output_file1, output_file2):
    """Validate paired-end FASTQ files."""
    # Check if first reads from both files have the same name
    with gzip.open(input_file1, 'rt') as fin1, gzip.open(input_file2, 'rt') as fin2:
        first_line1 = fin1.readline().strip()
        first_line2 = fin2.readline().strip()

        r1 = normalize_read_name(re.split(r" ", first_line1)[0])
        r2 = normalize_read_name(re.split(r" ", first_line2)[0])

        if r1 != r2:
            print("ERROR: First read names are different!", file=sys.stderr)
            print(f"       read 1: {r1}", file=sys.stderr)
            print(f"       read 2: {r2}", file=sys.stderr)
            sys.exit(1)

    # Count lines in both files
    n0 = count_gzip_lines(input_file1)
    n1 = count_gzip_lines(input_file2)

    # Validate line counts
    if n0 == n1:
        if n0 % 4 == 0:
            Path(output_file1).touch()
            Path(output_file2).touch()
        else:
            print(f"ERROR: Number of lines in FASTQ file is not a multiple of four: {n0}", file=sys.stderr)
            sys.exit(1)
    else:
        print("ERROR: Number of reads in mate files are different:", file=sys.stderr)
        print(f"File {input_file1} contains {n0} lines", file=sys.stderr)
        print(f"File {input_file2} contains {n1} lines", file=sys.stderr)
        sys.exit(1)


def main():
    """Main function to parse arguments and validate FASTQ files."""
    if len(sys.argv) < 4:
        print("Usage (Single-end): python validate_fastq.py single <input_file> <output_file>")
        print("Usage (Paired-end): python validate_fastq.py paired <input_file1> <input_file2> <output_file1> <output_file2>")
        sys.exit(1)

    mode = sys.argv[1]
    if mode == "single":
        if len(sys.argv) != 4:
            print("ERROR: Invalid arguments for single-end validation.")
            sys.exit(1)
        input_file = sys.argv[2]
        output_file = sys.argv[3]
        validate_fastq_single(input_file, output_file)

    elif mode == "paired":
        if len(sys.argv) != 6:
            print("ERROR: Invalid arguments for paired-end validation.")
            sys.exit(1)
        input_file1 = sys.argv[2]
        input_file2 = sys.argv[3]
        output_file1 = sys.argv[4]
        output_file2 = sys.argv[5]
        validate_fastq_paired(input_file1, input_file2, output_file1, output_file2)

    else:
        print("ERROR: Invalid mode. Use 'single' or 'paired'.")
        sys.exit(1)


if __name__ == "__main__":
    main()