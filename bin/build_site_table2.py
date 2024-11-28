from utils import read_fasta
from re import findall
import pandas as pd

args_msa = "/Users/stefanfrank/Projects/PanAdapt/results/split_reference_from_msa/ORF10_protein_0.pal2nal.no_ref"


def build_site_df(msa_file: str):
    site_df = pd.DataFrame()
    site_df.index.name = "Site"
    for identifier, sequence in read_fasta(msa_file):
        codons = findall("...", sequence)
        site_df[identifier] = codons
    site_df.insert(0, "Consensus_codon", site_df.mode(axis=1).iloc[:, 0])
    return site_df


def main():
    site_df = build_site_df(args_msa)
    site_df.to_csv("test.csv")


if __name__ == "__main__":
    main()
