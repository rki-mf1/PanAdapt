#!/usr/bin/env python
import argparse
from utils import read_fasta, write_fasta


def format_annotation(annotation, fasta):
    identifier, sequence = next(read_fasta(fasta))
    with open(annotation) as infile:
        lines = infile.readlines()
    lines.insert(1, f"##sequence-region {identifier.split()[0]} 1 {len(sequence)}\n")
    with open(annotation, "w") as outfile:
        outfile.writelines(lines)
        print("##FASTA", file=outfile)
        write_fasta(identifier, sequence, outfile)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--annotation", required=True)
    parser.add_argument("-f", "--fasta", required=True)
    args = parser.parse_args()
    format_annotation(args.annotation, args.fasta)


if __name__ == "__main__":
    main()
