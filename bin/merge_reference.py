#!/usr/bin/env python
import argparse
import fasta_utils
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-r", "--reference", required=True)
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep='\t', index_col=0)
    reference = fasta_utils.read_fasta(args.reference)
    reference_codons = reference[0].to_kmers(3)
    reference[0].translate()
    reference_aa = list(reference[0].protein)
    reference_position = -1
    reference_positions: list[int] = []
    for codon in reference_codons:
        if codon != "---":
            reference_position += 1
        reference_positions.append(reference_position)
    
    df.insert(loc=2, column="Reference_aa", value=reference_aa)
    df.insert(loc=1, column="Reference_codon", value=reference_codons)
    df.insert(loc=0, column="Reference_positition", value=reference_positions)

    df.to_csv(args.output, sep="\t")

if __name__ == "__main__":
    main()