#!/usr/bin/env python
import argparse
from utils import read_fasta, write_fasta

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--outfile")
    return parser.parse_args()

args = parse_args()
with open(args.outfile, "w") as outfile:
    for identifier, sequence in read_fasta(args.infile):
        if sequence[-3:] in ["TAA", "TAG", "TGA"]:
            sequence = sequence[:-3]
        write_fasta(identifier, sequence, outfile)