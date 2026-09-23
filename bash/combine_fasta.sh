#!/bin/bash

# Combine all mRNA FASTA files
for f in *mRNA*.fasta; do
    awk -v file="$f" '/^>/{print ">" file; next} 1' "$f"
done > all_species_combined_mRNA.fasta

# Combine all CDS FASTA files
for f in *CDS*.fasta; do
    awk -v file="$f" '/^>/{print ">" file; next} 1' "$f"
done > all_species_combined_CDS.fasta

# Combine all DNA FASTA files
for f in *DNA*.fasta; do
    awk -v file="$f" '/^>/{print ">" file; next} 1' "$f"
done > all_species_combined_DNA.fasta