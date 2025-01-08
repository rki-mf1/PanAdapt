#!/usr/bin/env python
import argparse
import fasta_utils
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-f", "--fubar", required=True)
    parser.add_argument("-r", "--reference", default=None)
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()

    fasta = fasta_utils.read_fasta(args.input)
    df = fasta.to_dataframe()
    df_aa = fasta.to_dataframe(k=1, mode="protein")
    fubar_df = pd.read_csv(args.fubar, sep="\t", index_col=0)
    shannon_entropy_codon = fasta.shannon_entropy(3)
    shannon_entropy_aa = fasta.shannon_entropy(k=1, mode="protein")
    fasta.get_consensus()
    consensus_codons = fasta.consensus.to_kmers(3)
    fasta.consensus.translate()
    consensus_aa = list(fasta.consensus.protein)

    if args.reference:
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

    fubar_df.insert(loc=0, column="Shannon_entropy_aa", value=shannon_entropy_aa)
    fubar_df.insert(loc=0, column="Shannon_entropy_codon", value=shannon_entropy_codon)
    fubar_df.insert(loc=0, column="Consensus_aa", value=consensus_aa)
    # fubar_df.insert(loc=0, column="Reference_aa", value=reference_aa)
    fubar_df.insert(loc=0, column="Consensus_codon", value=consensus_codons)
    # fubar_df.insert(loc=0, column="Reference_codon", value=reference_codons)
    # fubar_df.insert(loc=0, column="Reference_positition", value=reference_positions)

    df = fubar_df.merge(df, left_index=True, right_index=True, how="left")
    df = df.merge(
        df_aa, left_index=True, right_index=True, how="left", suffixes=("_codon", "_aa")
    )

    df.index.name = "Position"
    df.to_csv(args.output, sep="\t")


if __name__ == "__main__":
    main()
