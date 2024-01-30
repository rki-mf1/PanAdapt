#!/bin/bash
#SBATCH --job_name=PanAdapt
#SBATCH --time=5-00:00:00
#SBATCH --cpus-per-task=100
#SBATCH --mem-per-cpu=2G

WORKFLOW=$1
CONFIG=$2

nextflow -C ${CONFIG} run ${WORKFLOW} -profile slurm,conda