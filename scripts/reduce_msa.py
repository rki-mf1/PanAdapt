import argparse
from utils import read_fasta, write_fasta

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-t", "--tree")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-r", "--subtree")
    return parser.parse_args()

def extract_taxlabels(filename):
    taxlabels = []
    with open(filename, 'r') as file:
        lines = file.readlines()
        newick_tree = extract_newick_tree(lines)
        start, end = False, False
        for line in lines:
            line = line.strip()
            if line.startswith('TAXLABELS'):
                start = True
                continue
            if start and line == ";":
                end = True
            if start and not end:
                taxlabels.append(line.replace("'", ""))
    return taxlabels, newick_tree

def extract_newick_tree(lines):
    for line in lines:
        line = line.strip()
        if line.startswith('TREE'):
            return line.split('=')[1].strip()
    return None



def main():
    args = parse_args()
    taxlabels, newick_tree = extract_taxlabels(args.tree)
    with open(args.outfile, 'w') as outfile:
        for identifier, sequence in read_fasta(args.infile):
            if identifier in taxlabels:
                write_fasta(identifier, sequence, outfile)
    with open(args.subtree, 'w') as outfile:
        print(newick_tree, file=outfile)

if __name__ == '__main__':
    main()
