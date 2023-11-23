#!/usr/bin/env python
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--outfile")
    return parser.parse_args()


args = parse_args()

beb_line_found = False
amino_line_found = False
beb_table = []
with open(args.infile) as infile:
    for line in infile:
        if line.startswith("Bayes Empirical Bayes"):
            beb_line_found = True
            continue
        if beb_line_found:
            if line.startswith("(amino acids refer to 1st sequence:"):
                amino_line_found = True
                continue
        if beb_line_found and amino_line_found:
            if line.startswith("Positively selected sites"):
                break
            elif line != "\n":
                line = line[::-1]
                for char in ["+", "-", "(", ")"]:
                    line = line.replace(char, "", 1)
                line = line[::-1]
                line = line.split()
                beb_table.append(line)

column_names = (
    ["Position", "Amino_Acid"]
    + [f"Class_{i+1}_Prob" for i in range(11)]
    + ["Expected_Class", "Postmean_W", "Postmean_W_Error"]
)

with open(args.outfile, "w") as outfile:
    print("\t".join(column_names), file=outfile)
    for row in beb_table:
        print("\t".join(row), file=outfile)