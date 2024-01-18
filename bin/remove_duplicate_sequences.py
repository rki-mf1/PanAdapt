#!/usr/bin/env python
import argparse
import hashlib
import json
from utils import read_fasta, write_fasta


def remove_duplicate_sequences(infile, outfile_fasta, outfile_json, reference_id):
    hash_dict = {}
    ref_sequence_written = False

    with open(outfile_fasta, "w") as outfile:
        for identifier, sequence in read_fasta(infile):
            if reference_id in identifier:
                if not ref_sequence_written:
                    write_fasta(identifier, sequence, outfile)
                    ref_sequence_written = True
                continue

            hash_digest = hashlib.sha256(sequence.encode()).hexdigest()
            if hash_digest not in hash_dict:
                hash_dict[hash_digest] = [identifier]
                write_fasta(identifier, sequence, outfile)
            else:
                hash_dict[hash_digest].append(identifier)

    duplicates_dict = {
        idents[0]: idents[1:] for idents in hash_dict.values() if len(idents) > 1
    }

    with open(outfile_json, "w") as outfile:
        json.dump(duplicates_dict, outfile)


def main():
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
    parser.add_argument(
        "-r",
        "--reference",
        required=True,
        help="Part of the identifier of the reference sequence.",
    )
    args = parser.parse_args()
    remove_duplicate_sequences(args.infile, args.outfile, args.json, args.reference)


if __name__ == "__main__":
    main()
