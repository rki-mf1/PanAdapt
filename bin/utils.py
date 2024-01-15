#!/usr/bin/env python
import random
import argparse


def read_fasta(infile):
    with open(infile, "r") as file:
        identifier = None
        sequence = []
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if identifier is not None:
                    yield identifier, "".join(sequence)
                identifier = line[1:]
                sequence = []
            else:
                sequence.append(line)

        if identifier is not None:
            yield identifier, "".join(sequence)


def write_fasta(identifier, sequence, outfile, line_length=60):
    print(f">{identifier}", file=outfile)
    formatted_sequence = "\n".join(
        [sequence[i : i + line_length] for i in range(0, len(sequence), line_length)]
    )
    print(formatted_sequence, file=outfile)


def get_random_sequences(infile, n):
    reservoir = []
    for i, (identifier, sequence) in enumerate(read_fasta(infile)):
        if i < n:
            reservoir.append((identifier, sequence))
        else:
            j = random.randint(0, i)
            if j < n:
                reservoir[j] = (identifier, sequence)

    print(f"Processed {i + 1} sequences in total.")
    return reservoir


def main():
    parser = argparse.ArgumentParser(
        description="Utility script for handling FASTA files."
    )
    subparsers = parser.add_subparsers(dest="command", help="Commands")

    # Subparser for get_random_sequences
    parser_get_random = subparsers.add_parser(
        "get_random_sequences", help="Get random sequences from a FASTA file."
    )
    parser_get_random.add_argument(
        "-i", "--infile", required=True, help="Path to the input FASTA file."
    )
    parser_get_random.add_argument(
        "-o", "--outfile", required=True, help="Path to the output FASTA file."
    )
    parser_get_random.add_argument(
        "-n",
        "--number",
        type=int,
        required=True,
        help="Number of sequences to extract.",
    )

    args = parser.parse_args()

    if args.command == "get_random_sequences":
        random_sequences = get_random_sequences(args.infile, args.number)
        with open(args.outfile, "w") as outfile:
            for identifier, sequence in random_sequences:
                write_fasta(identifier, sequence, outfile)


if __name__ == "__main__":
    main()
