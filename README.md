# PanAdapt

## Introduction
PanAdapt is a `Nextflow` Pipeline designed to perform positive selection analysis on bacterial genomes. PanAdapt accepts genomes in fasta format as input and performs gene-level positive selection analysis using `BUSTED` and site-level positive selection analysis using `FUBAR`. 

## Installation
PanAdapt is available via github
```
git clone https://github.com/rki-mf1/PanAdapt
```
PanAdapt depends on `conda` for environment management. Please ensure that a version of conda is available on your system. PanAdapt was extensively tested with `miniconda`. The miniconda installation guide can be found [here](https://docs.conda.io/projects/miniconda/en/latest/). Using conda install and activate the `PanAdapt` environment with the following commands:
```
conda env create -f environment.yml
conda activate PanAdapt
```
This environment contains the packages `nextflow` and `bakta`. Now download a collection of databases for genome annotation with `bakta` using the following command:
```
bakta_db download --output resources --type light
```

## Running PanAdapt
Run the pipeline locally using the following command:
```
nextflow run main.nf --input_genomes path/to/genomes.fasta -profile local,conda 
```
Or run the pipeline via `sbatch` with the following command:
```
nextflow run main.nf --input_genomes path/to/genomes.fasta -profile slurm,conda 
```

## Output Files
By default all output are saved to the `results` directory in the base directory. Final results for the gene-level analysis can be found in the directory `extract_ppanggolon_results` and final results for site-level analysis are stored in the `calculate_shannon_entropy` directory. 