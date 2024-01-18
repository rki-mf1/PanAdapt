#!/usr/bin/env python
import argparse
import pandas as pd
import math


def calculate_shannon_entropy(row):
    unique_codons = row.value_counts()
    total = len(row)
    proportions = unique_codons / total
    entropy = -sum(p * math.log2(p) for p in proportions if p > 0)
    return entropy


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--outfile")
    args = parser.parse_args()
    df = pd.read_csv(args.infile, sep="\t")
    seq_columns = [col for col in df.columns if col.startswith("Seq_")]
    df["shannon_entropy"] = df[seq_columns].apply(calculate_shannon_entropy, axis=1)

    codon_to_aa = {
        "ATA": "I",
        "ATC": "I",
        "ATT": "I",
        "ATG": "M",
        "ACA": "T",
        "ACC": "T",
        "ACG": "T",
        "ACT": "T",
        "AAC": "N",
        "AAT": "N",
        "AAA": "K",
        "AAG": "K",
        "AGC": "S",
        "AGT": "S",
        "AGA": "R",
        "AGG": "R",
        "CTA": "L",
        "CTC": "L",
        "CTG": "L",
        "CTT": "L",
        "CCA": "P",
        "CCC": "P",
        "CCG": "P",
        "CCT": "P",
        "CAC": "H",
        "CAT": "H",
        "CAA": "Q",
        "CAG": "Q",
        "CGA": "R",
        "CGC": "R",
        "CGG": "R",
        "CGT": "R",
        "GTA": "V",
        "GTC": "V",
        "GTG": "V",
        "GTT": "V",
        "GCA": "A",
        "GCC": "A",
        "GCG": "A",
        "GCT": "A",
        "GAC": "D",
        "GAT": "D",
        "GAA": "E",
        "GAG": "E",
        "GGA": "G",
        "GGC": "G",
        "GGG": "G",
        "GGT": "G",
        "TCA": "S",
        "TCC": "S",
        "TCG": "S",
        "TCT": "S",
        "TTC": "F",
        "TTT": "F",
        "TTA": "L",
        "TTG": "L",
        "TAC": "Y",
        "TAT": "Y",
        "TAA": "_",
        "TAG": "_",
        "TGC": "C",
        "TGT": "C",
        "TGA": "_",
        "TGG": "W",
    }

    df_aa = df.replace(codon_to_aa)
    df["shannon_entropy_aa"] = df_aa[seq_columns].apply(
        calculate_shannon_entropy, axis=1
    )
    df.to_csv(args.outfile, sep="\t", index=False)


if __name__ == "__main__":
    main()
