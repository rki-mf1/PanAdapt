import argparse
import pandas as pd

def extract_positively_selected_sites(input_file, output_file):
    # Read the input file using pandas
    df = pd.read_csv(input_file, sep="\t")

    # Filter rows where 'Expected_Class' is 11
    filtered_df = df[df['Expected_Class'] == 11]

    # Write the filtered data to the output file
    filtered_df.to_csv(output_file, sep="\t", index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract positively selected sites from a given table.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input per site table.")
    parser.add_argument("-o", "--output", required=True, help="Path to the output file.")
    args = parser.parse_args()

    extract_positively_selected_sites(args.input, args.output)
