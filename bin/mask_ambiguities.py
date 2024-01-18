#!/usr/bin/env python
from utils import read_fasta, write_fasta
import re
import argparse


def transpose(infile):
    msa = read_fasta(infile)
    codons_by_seq = []
    for indentifier, seq in msa:
        codons_by_seq.append(re.findall(r"...", seq))
    codons_by_site = list(map(list, zip(*codons_by_seq)))
    transpose_dict = {}
    for i, site in enumerate(codons_by_site):
        count = {}
        for item in site:
            if all(n in "ACGT-" for n in item):
                if item in count:
                    count[item] += 1
                else:
                    count[item] = 1
        transpose_dict[i] = list(count.items())

    for key, value in transpose_dict.items():
        if len(value) == 1 and value[0][0] == "---":
            transpose_dict[key] = [("", value[0][1])]
    return transpose_dict


def consensus(transpose_dict):
    transpose_dict = {
        key: sorted(value, key=lambda x: x[1], reverse=True)
        for key, value in transpose_dict.items()
    }
    return [val[0][0] if val else "NNN" for val in transpose_dict.values()]


def replace_with_consensus(msa_file, consensus_seq, out_msa):
    msa = read_fasta(msa_file)
    with open(out_msa, "w") as outfile:  # Changed to 'w' mode
        for identifier, seq in msa:
            codons = re.findall(r"...", seq)
            for i, codon in enumerate(codons):
                if not all(n in "ACGT" for n in codon):
                    codons[i] = consensus_seq[i]  # Corrected this line
            updated_seq = "".join(codons)  # Join the codons back into a string
            write_fasta(identifier, updated_seq, outfile)  # Write the updated sequence


def main():
    parser = argparse.ArgumentParser(
        description="Process MSA and replace ambiguous codons with consensus."
    )
    parser.add_argument("-i", "--input", required=True, help="Input MSA file.")
    parser.add_argument("-o", "--output", required=True, help="Output MSA file.")
    args = parser.parse_args()

    # Compute consensus
    transposed_msa = transpose(args.input)
    cons_seq = consensus(transposed_msa)

    # Replace ambiguous codons in MSA with consensus codons
    replace_with_consensus(args.input, cons_seq, args.output)


if __name__ == "__main__":
    main()
