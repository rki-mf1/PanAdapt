import argparse
import pandas as pd
from utils import read_fasta, write_fasta

def build_lookup_table(csv_file):
    df = pd.read_csv(csv_file)

    # Identify the columns after "Avg group size nuc"
    identifier_cols = df.columns.tolist().index("Avg group size nuc") + 1

    # Create the dictionary
    lookup = {}

    for _, row in df.iterrows():
        gene = row['Gene']
        for identifier in row[identifier_cols:]:
            if pd.notna(identifier):
                lookup[identifier] = gene
    print(lookup)
    return lookup

def process_fasta_with_tsv(fasta_file, tsv_file, ref_id, outdir):
    lookup = build_lookup_table(tsv_file)
    for identifier, sequence in read_fasta(fasta_file):
        identifier = identifier.split()[0]
        if identifier in lookup:
            output_filename = lookup[identifier]
            if ref_id in identifier:
                with open(f"{outdir}/{output_filename}.ref", "a") as outfile:
                    write_fasta(identifier, sequence, outfile)
            else:
                with open(f"{outdir}/{output_filename}.fasta", "a") as outfile_1, open(f"{outdir}/{output_filename}.ref", "a") as outfile_2:
                    write_fasta(identifier, sequence, outfile_1)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--matrix")
    parser.add_argument("-s", "--split_cds")
    parser.add_argument("-r", "--ref_id")
    parser.add_argument("-o", "--outdir")
    args = parser.parse_args()
    process_fasta_with_tsv(args.split_cds, args.matrix, args.ref_id, args.outdir)

if __name__ == "__main__":
    main()
