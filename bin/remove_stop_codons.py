#!/usr/bin/env python
import argparse
from utils import read_fasta, write_fasta

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", required=True)
    parser.add_argument("-o", "--outfile", required=True)
    args = parser.parse_args()

    with open(args.outfile, "w") as outfile:
        for identifier, sequence in read_fasta(args.infile):
            if sequence.endswith(("TAA", "TAG", "TGA")):
                sequence = sequence[:-3]
            write_fasta(identifier, sequence, outfile)

if __name__ == "__main__":
    main()
