import argparse
import pandas as pd

def join_tables(per_site_file, beb_file, output_file):
    # Read the tables
    per_site_df = pd.read_csv(per_site_file, sep='\t', index_col=0)
    beb_df = pd.read_csv(beb_file, sep='\t', index_col=0)

    # Select only the desired columns from the beb table
    beb_df = beb_df[['Postmean_W', 'Postmean_W_Error']]

    # Merge the tables on the index (position)
    merged_df = per_site_df.join(beb_df, how='left')

    # Write the merged table to the output file
    merged_df.to_csv(output_file, sep='\t')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Join a per-site table with selected columns from a BEB table using pandas.")
    parser.add_argument('-i', '--input', required=True, help="Input per-site table file.")
    parser.add_argument('-b', '--beb', required=True, help="Input BEB table file.")
    parser.add_argument('-o', '--output', required=True, help="Output joined table file.")
    args = parser.parse_args()

    join_tables(args.input, args.beb, args.output)
