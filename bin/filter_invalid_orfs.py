#!/usr/bin/env python
import argparse

def remove_nonsense_sequences(input_file, output_file):
    with open(input_file) as infile, open(output_file, 'w') as outfile:
        valid_orf = True
        in_fasta_section = False

        for line in infile:
            if line.startswith('##FASTA'):
                in_fasta_section = True

            if in_fasta_section:
                outfile.write(line)
                continue

            if 'valid_ORF=True' in line:
                valid_orf = True
            elif 'valid_ORF=False' in line:
                valid_orf = False

            if valid_orf and 'ID=cds' in line:
                split_line = line.split('\t')
                complete_frame = (int(split_line[4]) - int(split_line[3]) + 1) % 3 == 0
                if complete_frame:
                    outfile.write(line)
            elif line.startswith('#'):
                outfile.write(line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()
    remove_nonsense_sequences(args.input, args.output)

