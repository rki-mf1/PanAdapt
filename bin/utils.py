#!/usr/bin/env python
import random
import argparse
import re


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


def translate(sequence: str, omit_terminal_stop: bool = False) -> str:
    translation_table = {
        "ATA": "I",
        "ATC": "I",
        "ATT": "I",
        "ATG": "M",
        "ACA": "T",
        "ACC": "T",
        "ACG": "T",
        "ACT": "T",
        "AAC": "N",
        "AAT": "N",
        "AAA": "K",
        "AAG": "K",
        "AGC": "S",
        "AGT": "S",
        "AGA": "R",
        "AGG": "R",
        "CTA": "L",
        "CTC": "L",
        "CTG": "L",
        "CTT": "L",
        "CCA": "P",
        "CCC": "P",
        "CCG": "P",
        "CCT": "P",
        "CAC": "H",
        "CAT": "H",
        "CAA": "Q",
        "CAG": "Q",
        "CGA": "R",
        "CGC": "R",
        "CGG": "R",
        "CGT": "R",
        "GTA": "V",
        "GTC": "V",
        "GTG": "V",
        "GTT": "V",
        "GCA": "A",
        "GCC": "A",
        "GCG": "A",
        "GCT": "A",
        "GAC": "D",
        "GAT": "D",
        "GAA": "E",
        "GAG": "E",
        "GGA": "G",
        "GGC": "G",
        "GGG": "G",
        "GGT": "G",
        "TCA": "S",
        "TCC": "S",
        "TCG": "S",
        "TCT": "S",
        "TTC": "F",
        "TTT": "F",
        "TTA": "L",
        "TTG": "L",
        "TAC": "Y",
        "TAT": "Y",
        "TAA": "*",
        "TAG": "*",
        "TGC": "C",
        "TGT": "C",
        "TGA": "*",
        "TGG": "W",
    }

    codons = re.findall("...", sequence)
    protein = "".join([translation_table.get(codon, "X") for codon in codons])
    if omit_terminal_stop and protein.endswith("*"):
        return protein[:-1]
    return protein


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
