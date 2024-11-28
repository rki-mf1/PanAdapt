from __future__ import annotations
from random import choices, randint
from re import findall
from textwrap import wrap
from itertools import groupby
from os.path import isdir
from statistics import multimode


class Sequence:

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

    def __init__(
        self, identifier: str = "", nucleotide: str = "", protein: str = ""
    ) -> None:
        self.identifier = identifier
        self.nucleotide = nucleotide
        self.protein = protein

    def __len__(self):
        return len(self.nucleotide)

    def to_string(self, mode: str = "nucleotide") -> str:
        if mode == "protein":
            wrapped_sequence = "\n".join(wrap(self.protein, 60))
            return f">{self.identifier}\n{wrapped_sequence}"
        elif mode == "nucleotide":
            wrapped_sequence = "\n".join(wrap(self.nucleotide, 60))
            return f">{self.identifier}\n{wrapped_sequence}"
        else:
            raise ValueError(
                f"Invalid mode {mode}. Mode must be either 'nucleotide' or 'protein'."
            )

    def random(self, length: int) -> None:
        self.nucleotide = "".join(choices("ACGT", k=length))

    def to_kmers(self, k: int = 3) -> list[str]:
        return findall("." * k, self.nucleotide)

    def translate(self, drop_terminal_stop: bool = False) -> None:
        codons = self.to_kmers()
        self.protein = "".join(
            [self.translation_table.get(codon, "X") for codon in codons]
        )
        if drop_terminal_stop:
            self.protein = self.protein.rstrip("*")

    def ambig_ratio(self) -> float:
        nucs = {"A", "T", "C", "G"}
        non_atcg_count = sum(1 for char in self.nucleotide if char not in nucs)
        return non_atcg_count / len(self)

    def contains_internal_stop(self) -> bool:
        return "*" in self.protein.rstrip("*")


class Fasta:

    def __init__(self, count: int = 0) -> None:
        self.records = [Sequence(identifier=f"Sequence_{i}") for i in range(count)]
        self.consensus = ""

    def __iter__(self):
        return iter(self.records)

    def __len__(self):
        return len(self.records)

    def __getitem__(self, index: int) -> Sequence:
        return self.records[index]

    def __setitem__(self, index: int, value: Sequence) -> None:
        self.records[index] = value

    def add_record(
        self, identifier: str, sequence: str, mode: str = "nucleotide"
    ) -> None:
        if mode == "protein":
            self.records.append(Sequence(identifier=identifier, protein=sequence))
        elif mode == "nucleotide":
            self.records.append(Sequence(identifier=identifier, nucleotide=sequence))

    def to_file(self, file_path: str, mode: str = "nucleotide") -> None:
        with open(file_path, "w") as outfile:
            for record in self.records:
                print(record.to_string(mode), file=outfile)

    def to_files(self, dir_path: str, mode: str = "nucleotide") -> None:
        if not isdir(dir_path):
            raise FileNotFoundError(f"Directory {dir_path} does not exist.")
        for record in self.records:
            with open(f"{dir_path}/{record.identifier}.fasta", "w") as outfile:
                print(record.to_string(mode), file=outfile)

    def to_string(self, mode: str = "nucleotide") -> str:
        return "\n".join(record.to_string(mode) for record in self.records)

    def random(self, length: int) -> None:
        for record in self.records:
            record.random(length)

    def translate(self, drop_terminal_stop: bool = False) -> None:
        for record in self.records:
            record.translate(drop_terminal_stop)

    def remove_duplicate_sequences(self) -> None:
        unique_records: dict[str, Sequence] = {}
        for record in self.records:
            unique_records[record.nucleotide] = record
        self.records = list(unique_records.values())

    def transpose(self, k: int) -> list[list[str]]:
        kmers_list = [record.to_kmers(k) for record in self.records]
        return list(zip(*kmers_list))

    def to_tsv(self, file_path: str, k: int) -> None:
        with open(file_path, "w") as outfile:
            for record in self.records:
                codons = record.to_kmers(k)
                codons.insert(0, record.identifier)
                print("\t".join(codons), file=outfile)

    def get_consensus(self) -> None:
        self.transposed_single = self.transpose(1)
        consensus_list: list[str] = []
        for nuc_list in self.transposed_single:
            valid_nucs = [nuc for nuc in nuc_list if nuc in "ACGT-"]
            consensus_list.append(multimode(valid_nucs)[0])
        self.consensus = "".join(consensus_list)

    def mask_ambiguities(self) -> None:
        if not self.consensus:
            self.get_consensus()
        for record in self.records:
            mask = (
                n2 if n1 not in "ACGT-" else n1
                for n1, n2 in zip(record.nucleotide, self.consensus)
            )
            record.nucleotide = "".join(mask)

    def get_random_sample(self, n: int) -> Fasta:
        sample_fasta = Fasta()
        for i, record in enumerate(self.records):
            if i < n:
                sample_fasta.add_record(record.identifier, record.nucleotide)
            else:
                j = randint(0, i)
                if j < n:
                    sample_fasta[j] = record
            print(i)
        return sample_fasta


def read_fasta(file_path: str) -> Fasta:
    fasta = Fasta()
    with open(file_path) as infile:
        lines = (line.strip() for line in infile)
        identifier = None
        for is_header, group in groupby(lines, lambda line: line.startswith(">")):
            if is_header:
                identifier = next(group)[1:].strip()
            else:
                sequence = "".join(group)
                if identifier is None:
                    raise ValueError("Sequence data encountered before header.")
                fasta.add_record(identifier, sequence)
        return fasta


def read_fasta_groups(file_path: str, lookup: dict[str, str]) -> dict[str, Fasta]:
    group_names = set(lookup.values())
    fasta_dict = {name: Fasta() for name in group_names}
    with open(file_path) as infile:
        lines = (line.strip() for line in infile)
        identifier = None
        for is_header, group in groupby(lines, lambda line: line.startswith(">")):
            if is_header:
                identifier = next(group)[1:].strip()
            else:
                sequence = "".join(group)
                if identifier is None:
                    raise ValueError("Sequence data encountered before header.")
                group_name = lookup[identifier]
                fasta_dict[group_name].add_record(identifier, sequence)
    return fasta_dict


def read_fasta_sample(file_path: str, n: int) -> Fasta:
    sample_fasta = Fasta()
    with open(file_path) as infile:
        lines = (line.strip() for line in infile)
        identifier = None
        for i, (is_header, group) in enumerate(
            groupby(lines, lambda line: line.startswith(">"))
        ):
            if is_header:
                identifier = next(group)[1:].strip()
            else:
                sequence = "".join(group)
                if identifier is None:
                    raise ValueError("Sequence data encountered before header.")
                if i // 2 < n:
                    sample_fasta.add_record(identifier, sequence)
                else:
                    j = randint(0, i)
                    if j < n:
                        sample_fasta.records[j] = Sequence(
                            identifier=identifier, nucleotide=sequence
                        )
    return sample_fasta


def main():
    # fasta = read_fasta("complete.fasta")
    sample = read_fasta_sample("complete.fasta", 25)
    for record in sample:
        record.nucleotide = record.nucleotide.replace("-", "N")
    sample.to_files("sc2_samples")
    # fasta.to_tsv("transpose.tsv", 3)

    # fasta.mask_ambiguities()
    # fasta.to_file("masked.fasta")


if __name__ == "__main__":
    main()
