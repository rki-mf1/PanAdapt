#!/usr/bin/env python
import argparse
import csv
import fasta_utils


def build_gene_families(matrix_path: str, fasta_path: str, ambig_cutoff: float) -> None:
    with open(matrix_path, mode="r") as infile:
        csv_reader = csv.reader(infile, delimiter="\t")
        lookup = {}
        next(csv_reader)
        for row in csv_reader:
            flattened_keys = [
                key.strip('"') for string in row[14:] for key in string.split()
            ]

            for key in flattened_keys:
                lookup[key] = row[0]

        fasta_dict = fasta_utils.read_fasta_groups(fasta_path, lookup)
        ambig_list: list[tuple[str, int]] = []
        stop_list: list[tuple[str, int]] = []
        for name, fasta in fasta_dict.items():
            fasta.records = [
                record
                for record in fasta.records
                if record.ambig_ratio() < ambig_cutoff
            ]
            ambig_list.append((name, len(fasta)))
            fasta.translate(
                drop_terminal_stop=True,
            )
            fasta.records = [
                record
                for record in fasta.records
                if not record.contains_internal_stop()
            ]
            stop_list.append((name, len(fasta)))
            if len(fasta) >= 2:
                fasta.to_file(f"{name}.fna")
                fasta.to_file(f"{name}.faa", mode="protein")

            with open(
                "ambig_filter_counts.tsv",
                "w",
            ) as outfile:
                writer = csv.writer(outfile, delimiter="\t")
                writer.writerow(["Gene", "Ambig_filter"])
                for row in ambig_list:
                    writer.writerow(row)

            with open(
                "stop_filter_counts.tsv",
                "w",
            ) as outfile:
                writer = csv.writer(outfile, delimiter="\t")
                writer.writerow(["Gene", "Stop_filter"])
                for row in stop_list:
                    writer.writerow(row)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--matrix", required=True)
    parser.add_argument("-f", "--fasta", required=True)
    parser.add_argument("-a", "--ambig_cutoff", type=float, default=1.0)
    args = parser.parse_args()
    build_gene_families(args.matrix, args.fasta, args.ambig_cutoff)


if __name__ == "__main__":
    main()
