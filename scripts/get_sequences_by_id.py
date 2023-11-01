import pandas as pd
from utils import read_fasta, write_fasta

df = pd.read_csv("sc2_input/2021-08_2022-03_random_no.tsv", sep="\t")
id_list = df["SEQUENCE.ID"].tolist()

with open("sc2_not_vaccinated.fasta", "w") as out_file:
    for identifier, sequence in read_fasta(
        "/mnt/wissdaten_sc2/RKI_nCoV-Lage/6.Datenaustausch/prod/complete_export/2023-10-24/complete.fasta"
    ):
        if identifier in id_list:
            write_fasta(identifier, sequence, out_file)
