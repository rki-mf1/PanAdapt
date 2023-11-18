import argparse
from utils import *


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-s", "--size", type=int)
    parser.add_argument("-n", "--number", type=int)
    return parser.parse_args()

def sanitize_sequence(sequence):
    return ''.join([base if base in 'ACGT' else 'N' for base in sequence])

def main():
    args = parse_args()
    random_sequences = get_random_sequences(args.infile, args.size, args.number)
    
    with open(args.outfile, "w") as out_file:
        for identifier, sequence in random_sequences:
            sanitized_sequence = sanitize_sequence(sequence)
            write_fasta(identifier, sanitized_sequence, out_file)


if __name__ == "__main__":
    main()