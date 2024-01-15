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
    df.to_csv(args.outfile, sep="\t", index=False)


if __name__ == "__main__":
    main()
