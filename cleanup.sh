#!/bin/bash
shopt -s extglob  # Enable extended globbing
rm -rf .nextflow*
rm -rf results/!(download_bakta_db)
rm -rf work/!(conda)
