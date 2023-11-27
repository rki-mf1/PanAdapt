#!/usr/bin/env python
import json
import pandas as pd
import argparse

def main(input_file, json_file, output_file):
    with open(json_file, 'r') as file:
        data = json.load(file)

    per_site_df = pd.read_csv(input_file, sep='\t', index_col=0)
    fubar_df = pd.DataFrame(data['MLE']['content']['0']).iloc[:, :6]
    fubar_df.columns = [header[0] for header in data['MLE']['headers']]
    merged_df = per_site_df.join(fubar_df, how='left')
    merged_df.to_csv(output_file, sep='\t')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert JSON data to TSV")
    parser.add_argument('-i', '--input_file', required=True, help="Input per-site table file.")
    parser.add_argument("-j", "--json_file", required=True, help="Path to the input JSON file")
    parser.add_argument("-o", "--output_file", required=True, help="Path to the output TSV file")
    args = parser.parse_args()

    main(args.input_file, args.json_file, args.output_file)
