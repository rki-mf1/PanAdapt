#!/usr/bin/env python
import argparse
from utils import read_fasta
import textwrap


def msa_to_tsv(infile, outfile_tsv):
    with open(outfile_tsv, "w") as outfile:
        for identifier, sequence in read_fasta(infile):
            codons = textwrap.wrap(sequence, 3)
            codons.insert(0, identifier)
            line = "\t".join(codons)
            print(line, file=outfile)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--outfile")
    args = parser.parse_args()
    msa_to_tsv(args.infile, args.outfile)


if __name__ == "__main__":
    main()
