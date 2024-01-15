#!/usr/bin/env python
import argparse
from utils import read_fasta, write_fasta
import json


def restore_duplicate_sequences(infile_fasta, infile_json, outfile_fasta):
    with open(infile_json, "r") as infile:
        duplicate_dict = json.load(infile)

    with open(outfile_fasta, "w") as outfile:
        for identifier, sequence in read_fasta(infile_fasta):
            duplicate_identifiers = duplicate_dict.get(identifier, [])
            write_fasta(identifier, sequence, outfile)
            for id in duplicate_identifiers:
                write_fasta(id, sequence, outfile)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-j", "--json")
    parser.add_argument("-o", "--outfile")
    args = parser.parse_args()
    restore_duplicate_sequences(args.infile, args.json, args.outfile)


if __name__ == "__main__":
    main()
