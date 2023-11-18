#!/usr/bin/env python
import argparse
import csv

def process_csv(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile)
        writer = csv.writer(outfile, delimiter='\t', quotechar='', quoting=csv.QUOTE_NONE)

        for row in reader:
            writer.writerow(row[:5])

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input CSV file")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file")
    args = parser.parse_args()
    process_csv(args.input, args.output)
