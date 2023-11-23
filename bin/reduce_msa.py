#!/usr/bin/env python
import argparse
from utils import read_fasta, write_fasta
import re

def get_taxlabels_from_newick(newick_tree_file):
    with open(newick_tree_file, 'r') as file:
        newick_tree = file.read().strip()
        labels = re.findall(r'[\w-]+_filtered_CDS_\d+', newick_tree)
        return labels

def main():
    parser = argparse.ArgumentParser(description="Reduce MSA based on tax labels in Newick tree.")
    parser.add_argument("-i", "--infile", required=True, help="Input FASTA file.")
    parser.add_argument("-t", "--tree", required=True, help="Input Newick tree file.")
    parser.add_argument("-o", "--outfile", required=True, help="Output reduced FASTA file.")
    args = parser.parse_args()

    taxlabels = get_taxlabels_from_newick(args.tree)
    with open(args.outfile, 'w') as outfile:
        for identifier, sequence in read_fasta(args.infile):
            if identifier in taxlabels:
                write_fasta(identifier, sequence, outfile)

if __name__ == '__main__':
    main()
