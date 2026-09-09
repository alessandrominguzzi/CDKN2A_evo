#!/bin/bash

# Launch this script exclusively inside the folder 'p16_candidates_blastn' 
# entire path here: ~/CDKN2A_evo/raw_data/sequences/p16_candidates_blastn
# this ensures that only p16 candidates found by blastn are included 

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