import argparse
import pandas as pd
from utils import read_fasta, write_fasta

def build_lookup_table(csv_file):
    df = pd.read_csv(csv_file)

    # Filter rows that contain the word 'persistent'
    df = df[df['Non-unique Gene name'].str.contains('persistent', case=False, na=False)]

    # Identify the columns after "Avg group size nuc"
    identifier_cols = df.columns.tolist().index("Avg group size nuc") + 1

    # Create the dictionary
    lookup = {}

    for _, row in df.iterrows():
        gene = row['Gene']
        annotation = row['Annotation'].replace(' ', '_')  # Replace spaces with underscores
        for identifier in row[identifier_cols:]:
            if pd.notna(identifier):
                lookup[identifier] = (annotation, gene)
    print(lookup)
    return lookup

def process_fasta_with_tsv(fasta_file, tsv_file, ref_id, outdir):
    lookup = build_lookup_table(tsv_file)
    for identifier, sequence in read_fasta(fasta_file):
        identifier = identifier.split()[0]
        if identifier in lookup:
            annotation, gene = lookup[identifier]
            # Only save the file if annotation is 'surface_glycoprotein'
            if annotation == 'surface_glycoprotein':
                output_filename = lookup[identifier][0]  # Use annotation (with underscores) as output filename
                if ref_id in identifier:
                    with open(f"{outdir}/{output_filename}.ref", "a") as outfile:
                        write_fasta(identifier, sequence, outfile)
                else:
                    with open(f"{outdir}/{output_filename}.fasta", "a") as outfile:
                        write_fasta(identifier, sequence, outfile)

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
