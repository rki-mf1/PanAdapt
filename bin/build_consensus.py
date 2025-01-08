#!/usr/bin/env python
import argparse
import fasta_utils


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()

    fasta = fasta_utils.read_fasta(args.input)
    fasta.get_consensus()
    with open(args.output, "w") as outfile:
        print(fasta.consensus.to_string(), file=outfile)


if __name__ == "__main__":
    main()
