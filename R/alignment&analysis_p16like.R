# ==============================================================================
## ALIGNMENT OF p16 CANDIDATES ##
# ==============================================================================

# This script aims to align the mRNA and CDS sequences obtained from blastn 
# using the 'human_p16_CDS.fasta" query, in order to assess sequence   
# conservation and evaluate the plausibility of their homology. 

# ==============================================================================
# 1. mRNA Alignment
# ==============================================================================

# Loading libraries
library(Biostrings)
library(DECIPHER)
library(ggplot2)

# Setting paths

fasta_dir <- path.expand("~/CDKN2A_evo/raw_data/sequences/p16_candidates_blastn")

metadata_dir <- path.expand("~/CDKN2A_evo/raw_data/metadata/p16_blast_sheet.xlsx")

setwd(fasta_dir)

# Obtaining alignment

mrna_seqs <- readDNAStringSet("p16like_all_species_combined_mRNA.fasta")
aligned_mRNA <- AlignSeqs(mrna_seqs, verbose = FALSE)

#writeXStringSet(aligned_mRNA,  filepath = "/mnt/scratch/home/alessandro/CDKN2A_evo/processed_data/p16like_aligned_mRNA.fasta")

# ==============================================================================
# 2. CDS Alignment and Analysis
# ==============================================================================

## 2.1 Reading CDS sequences

cds_seqs <- readDNAStringSet("p16like_all_species_combined_CDS.fasta")
   

## 2.2 ORF integrity check

#to verify each CDS has valid ORFs: 
# 1) lenght x3
# 2) no premature stop codons

orf_length_ok <- (width(cds_seqs) %% 3) == 0

# translate()
translated <- AAStringSet(rep("", length(cds_seqs)))
translated[orf_length_ok] <- suppressWarnings(
  translate(cds_seqs[orf_length_ok], if.fuzzy.codon = "solve")
)

# Counting stop codons (*) in translated proteins, excluding the final
# stop codon
count_premature_stops <- function(aa_seq) {
  aa_chars <- strsplit(as.character(aa_seq), "")[[1]]
  # remove last character if stop codon (normal)
  if (length(aa_chars) > 0 && aa_chars[length(aa_chars)] == "*") {
    aa_chars <- aa_chars[-length(aa_chars)]
  }
  sum(aa_chars == "*")
}

premature_stops <- rep(NA_integer_, length(cds_seqs))
premature_stops[orf_length_ok] <- vapply(
  translated[orf_length_ok],
  count_premature_stops,
  integer(1)
)

orf_check <- data.frame(
  sequence_name        = names(cds_seqs),
  length_nt            = width(cds_seqs),
  multiple_of_3        = orf_length_ok,
  premature_stop_count = premature_stops,
  orf_ok               = orf_length_ok & (premature_stops == 0)
)

orf_check

## 2.3 Alignment 

aligned_cds <- AlignSeqs(cds_seqs, verbose = FALSE)
#writeXStringSet(aligned_cds, filepath = "/mnt/scratch/home/alessandro/CDKN2A_evo/processed_data/p16like_aligned_CDS.fasta")

## 2.4 QC and EDA

# Checking lengths
width(cds_seqs)
width(aligned_cds)


