import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge four TSV files based on the 'Gene' column."
    )
    parser.add_argument("file1", help="Path to the first TSV file.")
    parser.add_argument("file2", help="Path to the second TSV file.")
    parser.add_argument("file3", help="Path to the third TSV file.")
    # parser.add_argument("file4", help="Path to the fourth TSV file.")
    parser.add_argument(
        "-o", "--outfile", required=True, help="Path to the merged output TSV file."
    )
    return parser.parse_args()


def main(args):
    df1 = pd.read_csv(args.file1, sep="\t")
    df2 = pd.read_csv(args.file2, sep="\t")
    df3 = pd.read_csv(args.file3, sep="\t")
    # df4 = pd.read_csv(args.file4, sep="\t")
    merged_df = (
        df1.merge(df2, on="Gene", how="outer")
        .merge(df3, on="Gene", how="outer")
        # .merge(df4, on="Gene", how="outer")
    )
    merged_df = merged_df.astype(str)
    merged_df.to_csv(args.outfile, sep="\t", index=False)


if __name__ == "__main__":
    args = parse_args()
    main(args)
