import argparse

def remove_nonsense_sequences(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        write_cds = True
        for line in infile:
            if line.startswith('#') or 'ID=cds' in line:
                if write_cds:
                    outfile.write(line)
            elif 'ID=gene' in line:
                if 'valid_ORF=False' not in line:
                    write_cds = True
                else:
                    write_cds = False

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove nonsense sequences from GFF annotation.")
    parser.add_argument("-i", "--input", help="Input GFF annotation file")
    parser.add_argument("-o", "--output", help="Output GFF file after removing nonsense sequences")
    args = parser.parse_args()
    
    if args.input and args.output:
        remove_nonsense_sequences(args.input, args.output)
    else:
        print("Please provide both input and output file paths.")
