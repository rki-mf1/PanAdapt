#!/usr/bin/env python

import pandas as pd
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_matrix", required=True)
    parser.add_argument("-a", "--ambig_filter_counts", required=True)
    parser.add_argument("-s", "--stop_filter_counts", required=True)
    parser.add_argument("-d", "--merged_dups_filter", required=True)
    parser.add_argument("-f", "--merged_fubar_stats", required=True)
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()

    df_matrix = pd.read_csv(args.input_matrix, sep="\t")
    df_ambig = pd.read_csv(args.ambig_filter_counts, sep="\t")
    df_stop = pd.read_csv(args.stop_filter_counts, sep="\t")
    df_dups = pd.read_csv(args.merged_dups_filter, sep="\t")
    df_fubar = pd.read_csv(args.merged_fubar_stats, sep="\t")

    df_merged = df_matrix.merge(df_ambig, on="Gene", how="left")
    df_merged = df_merged.merge(df_stop, on="Gene", how="left")
    df_merged = df_merged.merge(df_dups, on="Gene", how="left")
    df_merged = df_merged.merge(df_fubar, on="Gene", how="left")

    df_merged["Dups_filter"] = df_merged["Dups_filter"].fillna(0)
    df_merged["Dups_filter"] = df_merged["Dups_filter"].astype(int)

    df_merged.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
