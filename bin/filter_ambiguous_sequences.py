#!/usr/bin/env python
import argparse
from utils import read_fasta, write_fasta


def calculate_amig_ratio(seq):
    non_atgc_count = sum(1 for char in seq if char not in "ATGC")
    total_count = len(seq)
    return non_atgc_count / total_count


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file")
    parser.add_argument("-o", "--output_file")
    args = parser.parse_args()
    with open(args.output_file, "w") as outfile:
        for identifier, sequence in read_fasta(args.input_file):
            if calculate_amig_ratio(sequence) < 0.1:
                write_fasta(identifier, sequence, outfile)


if __name__ == "__main__":
    main()
