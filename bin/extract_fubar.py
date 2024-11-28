#!/usr/bin/env python
import argparse
import json
import pandas as pd


def extract_fubar(fubar_json: str, gene: str, outfile: str, fubar_stats: str) -> None:
    with open(fubar_json, "r") as file:
        data = json.load(file)
    fubar_df = pd.DataFrame(data["MLE"]["content"]["0"]).iloc[:, :6]
    fubar_df.columns = [header[0] for header in data["MLE"]["headers"]]
    fubar_df.to_csv(outfile, sep="\t")

    min_value = fubar_df["Prob[alpha<beta]"].min()
    mean_value = fubar_df["Prob[alpha<beta]"].mean()
    max_value = fubar_df["Prob[alpha<beta]"].max()

    result_df = pd.DataFrame(
        {
            "Gene": [gene],
            "min_Prob[alpha<beta]": [min_value],
            "mean_Prob[alpha<beta]": [mean_value],
            "max_Prob[alpha<beta]": [max_value],
        }
    )

    result_df.to_csv(fubar_stats, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-g", "--gene", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-s", "--stats", required=True)
    args = parser.parse_args()
    extract_fubar(args.input, args.gene, args.output, args.stats)


if __name__ == "__main__":
    main()
