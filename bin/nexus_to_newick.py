#!/usr/bin/env python
import argparse

def extract_newick_tree(filename):
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('TREE'):
                tree = line.split('=')[1].strip()
                if not tree.startswith('(') and tree.endswith(';'):
                    tree = '(' + tree.rstrip(';') + '):'
                return tree
    return None

def main():
    parser = argparse.ArgumentParser(description="Convert NEXUS to Newick format.")
    parser.add_argument("-i", "--infile", required=True, help="Input NEXUS file.")
    parser.add_argument("-o", "--outfile", required=True, help="Output Newick file.")
    args = parser.parse_args()

    newick_tree = extract_newick_tree(args.infile)
    if newick_tree:
        with open(args.outfile, 'w') as outfile:
            print(newick_tree, file=outfile)

if __name__ == '__main__':
    main()
