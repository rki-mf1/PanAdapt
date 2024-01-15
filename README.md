# UpDownSelect

## seqkit
Splits the input multi-sequence fasta file into multiple single-sequence fasta files. The file names of the output files are equal to the individual sequence headers in the input file followed by a .fasta file extension.

## liftoff
Applies reference annotation to the single-sequence fasta files. '##FASTA' and single-sequence input fasta are appended to output .gff file.

## filter_invalid_orfs
Takes individual modified liftoff .gff file and removes all annotated features, that do not contain CDS with a length divisable by 3. The filtered file is saved with .ggf_filtered as the file extension.

## ppangolin
Takes all modified and filtered gff-files and builds a pangenome. Produces two relevant files, `all_genes.fna`, which contains every coding sequence from the input gff files in fasta format and `matrix.csv`, which combinesd the fasta headers from `all_genes.fna` into gene_families.

## extract_ppanggolin_results
Takes `matrix.csv` and extracts the columns `Gene`, `Non-unique`, `Gene name`, `Annotation`, `No. isolates`, `No. sequences`. The output is saved as a tsv-file.

## remove_stop_codons
Takes `all_genes.fna` as input. Every sequence ending in a stop codon has that codon removed. All truncated sequences are written to a single multi-sequence fasta file `all_genes.fasta`

## build_gene_families
Takes `matrix.csv` and extracts the fasta headers for each gene family. Builds one gene family fasta file for each row in the input file. The files are named according to the fasta header in the `Gene` column of `matrix.csv`

## remove_duplicate_sequences
Takes a gene family fasta-files and finds identical sequences. Only the first sequence of each duplicate group is written to a new output fasta file. Additional duplicate identifiers are collected in a json output file.

## filter_min_seq_count
Takes a gene family fasta file and counts the number of sequences contained within it. If the file contains more than 2 sequences, the gene famility is considered for further analysis. Files with 2 or less sequences are discarded.

## translate_gene_familise
Takes a filtered gene family fasta file and translated the containing sequences into amino acids using the standard translation table.

## mafft