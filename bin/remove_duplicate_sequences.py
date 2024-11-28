#!/usr/bin/env python
import fasta_utils
import argparse
import csv


def remove_duplicates_sequences(infile: str, outfile: str, gene: str):
    fasta = fasta_utils.read_fasta(infile)
    fasta.remove_duplicate_sequences()
    if len(fasta) >= 2:
        fasta.to_file(outfile)
    with open(f"{gene}_dups_filter_counts.tsv", "w") as outfile2:
        writer = csv.writer(outfile2, delimiter="\t")
        writer.writerow(["Gene", "Dups_filter"])
        writer.writerow([gene, len(fasta)])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", required=True)
    parser.add_argument("-o", "--outfile", required=True)
    parser.add_argument("-g", "--gene", required=True)
    args = parser.parse_args()
    remove_duplicates_sequences(args.infile, args.outfile, args.gene)
