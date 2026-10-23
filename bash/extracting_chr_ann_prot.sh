#!/bin/bash

# Chromosome accessions are stored in the metadata (Excel files).
# The following commands are only examples (human chromosome 9)
# for annotation extraction (gff file) and protein extraction (faa file).
# The same procedure must be performed for every species included in the study.

# Extract annotations for the chromosome of interest 
awk '$1 == "#" || $1 == "NC_000009.12"' human_genomic.gff > human_chr9.gff

# Extract proteins for the chromosome of interest 
seqkit grep -r -p "NC_000009.12" human.faa > human_ch9.faa