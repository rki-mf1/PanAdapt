#!/usr/bin/env python
import argparse
from utils import read_fasta, write_fasta

def split_reference_from_msa(args):
    with open(args.output, 'w') as outfile, open(args.ref_output, 'w') as ref_outfile:
        for identifier, sequence in read_fasta(args.input):
            if args.ref_id in identifier:
                write_fasta(identifier, sequence, ref_outfile)
            else:
                write_fasta(identifier, sequence, outfile)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split reference sequence from MSA based on ref_id")
    parser.add_argument("-i", "--input", required=True, help="Path to the input MSA file")
    parser.add_argument("-r", "--ref_id", required=True, help="Reference ID to search for in the MSA headers")
    parser.add_argument("-o", "--output", required=True, help="Path to the output MSA file without reference sequence")
    parser.add_argument("-u", "--ref_output", required=True, help="Path to the output file containing the reference sequence")

    args = parser.parse_args()
    split_reference_from_msa(args)
