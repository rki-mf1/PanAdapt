import argparse
import json
import re
from collections import Counter
from utils import read_fasta


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-d", "--identifier")
    return parser.parse_args()


def build_single_subs(codon):
    bases = "ACGT"
    stop_codons = ["TAA", "TAG", "TGA"]
    substitutions = []
    for i, base in enumerate(codon):
        for alt_base in bases:
            if alt_base != base:
                mutated_codon = list(codon)
                mutated_codon[i] = alt_base
                mutated_codon = "".join(mutated_codon)
                if mutated_codon not in stop_codons:
                    substitutions.append("".join(mutated_codon))

    return substitutions


def build_synonymy_dict(codon_table):
    with open(codon_table) as infile:
        codon_dict = json.load(infile)

    synonymy_dict = {}
    for key, value in codon_dict.items():
        subs_list = build_single_subs(key)
        syn_subs = [sub for sub in subs_list if codon_dict[sub] == codon_dict[key]]
        nonsyn_subs = [sub for sub in subs_list if sub not in syn_subs]
        synonymy_dict[key] = {"synonymous": syn_subs, "non-synonymous": nonsyn_subs}

    return synonymy_dict


def extract_site_patterns(infile):
    alignment = read_fasta(infile)
    codons_per_seq = []
    for id, seq in alignment:
        codons_per_seq.append(re.findall(r"...", seq))
    site_patterns = [site for site in zip(*codons_per_seq)]
    site_patterns = [dict(Counter(site)) for site in site_patterns]
    site_patterns = [
        dict(sorted(site.items(), key=lambda item: item[1], reverse=True))
        for site in site_patterns
    ]
    return site_patterns


def calculate_dNdS(infile, outfile, identifier, codon_table):
    synonymy_dict = build_synonymy_dict(codon_table)
    site_patterns = extract_site_patterns(infile)
    with open(codon_table) as file:
        codon_dict = json.load(file)
    syn_count = 0.001
    nonsyn_count = 0
    syn_sites = 0
    nonsyn_sites = 0
    for site in site_patterns:
        ref_codon = None
        for key, value in site.items():
            if key in codon_dict.keys():
                if ref_codon == None:
                    ref_codon = key
                    syn_sites += len(synonymy_dict[key]["synonymous"])
                    nonsyn_sites += len(synonymy_dict[key]["non-synonymous"])
                else:
                    syn_count += value * (codon_dict[key] == codon_dict[ref_codon])
                    nonsyn_count += value * (codon_dict[key] != codon_dict[ref_codon])
    dN = nonsyn_count / nonsyn_sites
    dS = syn_count / syn_sites
    with open(outfile, "a") as file:
        print(f"{identifier}\t{dN/dS}", file=file)


args = parse_args()
calculate_dNdS(
    args.infile,
    args.outfile,
    args.identifier,
    "resources/default_codon_table.json",
)
