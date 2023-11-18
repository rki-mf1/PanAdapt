#!/usr/bin/env python
import argparse
import pandas as pd
from utils import read_fasta, write_fasta

def build_lookup_table(csv_file):
    df = pd.read_csv(csv_file)
    return {row[col]: row['Annotation'].replace(' ', '_')
            for _, row in df.iterrows()
            for col in df.columns[df.columns.get_loc("Avg group size nuc") + 1:]}

def process_fasta_with_tsv(fasta_file, tsv_file, ref_id):
    lookup = build_lookup_table(tsv_file)
    for identifier, sequence in read_fasta(fasta_file):
        key = identifier.split()[0]
        if key in lookup:
            with open(f"{lookup[key]}.fasta", "a") as outfile:
                write_fasta(identifier, sequence, outfile)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--matrix", required=True)
    parser.add_argument("-s", "--split_cds", required=True)
    parser.add_argument("-r", "--ref_id", required=True)
    args = parser.parse_args()
    process_fasta_with_tsv(args.split_cds, args.matrix, args.ref_id)

if __name__ == "__main__":
    main()
