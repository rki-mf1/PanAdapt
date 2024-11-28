#!/usr/bin/env python
import argparse
import pandas as pd


def reformat_matrix(input_filepath, output_filepath):
    matrix = pd.read_csv(input_filepath)
    matrix = matrix.sort_values(
        by=["Annotation", "No. sequences"], ascending=[False, False]
    )
    matrix["Gene"] = (
        matrix["Annotation"] + "_" + matrix.groupby("Annotation").cumcount().astype(str)
    )
    matrix["Gene"] = matrix["Gene"].str.replace(r"[^A-Za-z0-9]", "_", regex=True)
    matrix.to_csv(output_filepath, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()
    reformat_matrix(args.input, args.output)


if __name__ == "__main__":
    main()
