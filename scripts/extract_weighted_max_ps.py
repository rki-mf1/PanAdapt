import argparse
import pandas as pd
import os


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--indir",
        required=True,
        help="Path to the input directory containing TSV files.",
    )
    parser.add_argument(
        "-o", "--outfile", required=True, help="Path to the output file."
    )
    return parser.parse_args()


def process_file(infile, outfile):
    df = pd.read_csv(infile, sep="\t")
    df["Weighted_PS"] = df["Postmean_W"] * df["Class_11_Prob"]
    base_name = os.path.basename(infile).split(".")[0]
    max_val = df["Weighted_PS"].max()

    # Format the max value to have 3 decimal places
    formatted_max_val = "{:.3f}".format(max_val)

    with open(outfile, "a") as f:
        f.write(f"{base_name}\t{formatted_max_val}\n")


def main(args):
    if not os.path.isdir(args.indir):
        print(f"Error: {args.indir} is not a directory.")
        return

    for filename in os.listdir(args.indir):
        if filename.endswith(".tsv"):
            infile = os.path.join(args.indir, filename)
            process_file(infile, args.outfile)


if __name__ == "__main__":
    args = parse_args()
    main(args)
