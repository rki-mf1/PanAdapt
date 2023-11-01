import argparse
import hashlib
import json
from utils import read_fasta, write_fasta


def parse_args():
    parser = argparse.ArgumentParser(
        description="Remove duplicate sequences from a fasta file."
    )
    parser.add_argument("-i", "--infile", required=True, help="Input fasta file.")
    parser.add_argument(
        "-o", "--outfile", required=True, help="Output fasta file without duplicates."
    )
    parser.add_argument(
        "-j",
        "--json",
        required=True,
        help="Output json file mapping original identifiers to duplicates.",
    )
    return parser.parse_args()


def remove_duplicate_sequences(infile, outfile_fasta, outfile_json):
    hash_dict = {}
    count_duplicates = 0
    total_sequences = 0

    with open(outfile_fasta, "w") as outfile:
        for identifier, sequence in read_fasta(infile):
            total_sequences += 1
            hash = hashlib.sha256(sequence.encode()).hexdigest()
            if hash not in hash_dict:
                write_fasta(identifier, sequence, outfile)
                hash_dict[hash] = [identifier]
            else:
                count_duplicates += 1
                hash_dict[hash].append(identifier)

    # Create a dictionary mapping the first identifier of a sequence to its duplicates
    hash_dict = {v[0]: v[1:] for _, v in hash_dict.items()}

    with open(outfile_json, "w") as outfile:
        json.dump(hash_dict, outfile)

    # Print messages
    print(f"Processing fasta file: {infile}")
    print(f"Total sequences processed: {total_sequences}")
    print(f"Total unique sequences: {total_sequences - count_duplicates}")
    print(f"Number of duplicate sequences: {count_duplicates}")
    print(f"Output fasta file: {outfile_fasta}")
    print(f"Output json file containing duplicates: {outfile_json}")
    print("Finished processing.")


if __name__ == "__main__":
    args = parse_args()
    remove_duplicate_sequences(args.infile, args.outfile, args.json)
