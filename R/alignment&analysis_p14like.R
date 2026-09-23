# ==============================================================================
## ALIGNMENT OF p14 CANDIDATES ##
# ==============================================================================

# This script aims to align the mRNA and CDS sequences obtained from blastn 
# using the 'human_p14_CDS.fasta" query, in order to assess sequence   
# conservation and evaluate the plausibility of their homology. 


# ==============================================================================
# 1. mRNA Alignment
# ==============================================================================

#Loading libraries

library(Biostrings)
library(DECIPHER)
library(ggplot2)
library(readxl)

#Setting paths and loading metadata

metadata_dir <- path.expand("~/CDKN2A_evo/raw_data/metadata/p14_blastn_sheet.xlsx")
metadata <- read_excel(metadata_dir)

fasta_dir <- path.expand("~/CDKN2A_evo/raw_data/sequences/p14_candidates_blastn")
setwd(fasta_dir)

# Loading metadata

metadata <- read_excel(metadata_dir)

#obtaining alignment

mrna_seqs <- readDNAStringSet("p14like_all_species_combined_mRNA.fasta")
aligned_mRNA <- AlignSeqs(mrna_seqs, verbose = FALSE)

#writeXStringSet(aligned_mRNA,  filepath = "/mnt/scratch/home/alessandro/CDKN2A_evo/processed_data/p14like_aligned_mRNA.fasta")


# ==============================================================================
# 2. CDS Alignment and Analysis
# ==============================================================================

## 2.1 Reading CDS sequences

cds_seqs <- readDNAStringSet("p14like_all_species_combined_CDS.fasta") 

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
#writeXStringSet(aligned_cds, filepath = "aligned_CDS.fasta")

## 2.4 QC and EDA

# Checking lengths
width(cds_seqs)
width(aligned_cds)

# Creating a matrix 
mat <- as.matrix(aligned_cds)
n_seqs <- nrow(mat)   

table(mat)

# Creating a species table per sequence
species_patterns <- c(
  "Human"      = "human",
  "Anolis"     = "anolis",
  "Blue_whale" = "bluewhale",   
  "Chicken"    = "chicken",
  "Mouse"      = "mouse",
  "Opossum"    = "opossum",
  "Pelodiscus" = "pelodiscus",
  "Platypus"   = "platypus",
  "Xenopus"    = "xenopus",
  "Zebrafish"  = "zebrafish"
)

species <- vapply(names(cds_seqs), function(nm) {
  hit <- names(species_patterns)[vapply(
    species_patterns,
    function(pat) grepl(pat, nm, ignore.case = TRUE),
    logical(1)
  )]
  if (length(hit) == 0) NA_character_ else hit[[1]]
}, character(1))

stopifnot(
  "Some sequences are not assigned to a species: control species_patterns / fasta header " =
    !any(is.na(species))
)


data.frame(sequence_name = names(cds_seqs), species = species)

species_table <- table(species)
species_table


# Evaluate contribution of same species sequences
multi_species <- names(species_table[species_table > 1])

for (sp in multi_species) {
  sp_mat <- as.matrix(aligned_cds[species == sp])
  modal_count <- apply(sp_mat, 2, function(x) {
    max(table(x))
  })
  
  cat("\n", sp, "\n")
  print(table(modal_count))
}

# Investigating gaps
gap_frequency <- apply(mat, 2, function(x) {
  mean(x == "-")
})

summary(gap_frequency)
table(gap_frequency)

gap_counts <- round(gap_frequency * n_seqs)  

gap_table <- data.frame(
  gap_count = as.integer(names(table(gap_counts))),
  nucleotide_count = n_seqs - as.integer(names(table(gap_counts))),  
  column_count = as.integer(table(gap_counts))
)

gap_table$gap_frequency <- gap_table$gap_count / n_seqs  
gap_table$column_percentage <- (gap_table$column_count / length(gap_frequency)) * 100
gap_table


# Descriptive analysis of the alignment
consensus <- apply(mat, 2, function(x) {
  nucleotides <- x[x %in% c("A", "C", "G", "T")]
  
  if (length(nucleotides) == 0) {
    return(NA_character_)
  }
  names(which.max(table(nucleotides)))
})

agreement <- apply(mat, 2, function(x) {
  nucleotides <- x[x %in% c("A", "C", "G", "T")]
  
  if (length(nucleotides) == 0) {
    return(NA_real_)
  }
  
  max(table(nucleotides)) / length(nucleotides)
})

occupancy <- 1 - gap_frequency

consensus[occupancy < 0.50] <- NA_character_ #0.5 occupancy threshold defined 

conservation_table <- data.frame(
  position = 1:length(consensus),
  consensus = consensus,
  agreement = agreement,
  occupancy = occupancy,
  gap_frequency = gap_frequency
)

#head(conservation_table, 20)

na_count <- sum(is.na(conservation_table$consensus)) 
na_percentage <- mean(is.na(conservation_table$consensus)) * 100

conservation_table$agreement_percent <- conservation_table$agreement * 100
conservation_table$occupancy_percent <- conservation_table$occupancy * 100

## 2.3 GRAPHS

#common theme
p14_theme <- theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

# Nucleotide conservation across the CDS alignment
ggplot(conservation_table, aes(x = position, y = agreement_percent)) +
  geom_line(linewidth = 0.5, color = "steelblue", alpha = 0.5) +
  geom_point(
    aes(color = occupancy_percent >= 50),
    size = 1.2,
    alpha = 0.5
  ) +
  geom_smooth(
    method = "loess",
    se = FALSE,
    color = "blue",
    linewidth = 1
  ) +
  geom_hline(
    yintercept = 50,
    linetype = "dashed",
    color = "grey40"
  ) +
  scale_color_manual(
    values = c("FALSE" = "tomato", "TRUE" = "darkgreen"),
    labels = c("FALSE" = "Occupancy < 50%", "TRUE" = "Occupancy ≥ 50%"),
    name = "Occupancy"
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20)
  ) +
  labs(
    title = "Nucleotide conservation across the CDS alignment",
    x = "Alignment position",
    y = "Agreement (%)"
  ) +
  p14_theme +
  theme(legend.position = "top")


#Sequence occupancy across the CDS alignment
ggplot(conservation_table, aes(x = position, y = occupancy_percent)) +
  geom_line(linewidth = 0.5, alpha = 0.5) +
  geom_smooth(
    method = "loess",
    se = FALSE,
    color = "black",
    linewidth = 1
  ) +
  geom_hline(
    yintercept = 50,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20)
  ) +
  labs(
    title = "Sequence occupancy across the CDS alignment",
    x = "Alignment position",
    y = "Occupancy (%)"
  ) +
  p14_theme

# Creating a conservation score: a measure that takes into account at the same
# time agreement and occupancy 
conservation_table$conservation_score <-
  conservation_table$agreement * conservation_table$occupancy * 100

#Conservation tracked based on the conservation score
ggplot(conservation_table, aes(
  x = position,
  y = 1,
  fill = conservation_score
)) +
  geom_tile(height = 1) +
  scale_fill_viridis_c(
    name = "Conservation\nscore (%)",
    limits = c(0, 100)
  ) +
  scale_y_continuous(
    breaks = NULL
  ) +
  labs(
    title = "Conservation track across the CDS alignment",
    x = "Alignment position",
    y = NULL
  ) +
  p14_theme +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

##
ggplot(
  conservation_table,
  aes(x = occupancy_percent, y = agreement_percent)
) +
  geom_point(alpha = 0.5, size = 2) +
  geom_vline(
    xintercept = 50,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20)
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20)
  ) +
  labs(
    title = "Agreement and sequence occupancy across CDS positions",
    x = "Occupancy (%)",
    y = "Agreement (%)"
  ) +
  p14_theme

#Heatmap species x  position
# --------------------------------------------------------------------------
# species per line to see where the gaps are, useful to investigate if gaps
# clustered in a few species or if they are randomly distributed along several
# ones. 


heatmap_long <- data.frame(
  sequence_name = rep(names(cds_seqs), times = ncol(mat)),
  species       = rep(species, times = ncol(mat)),
  position      = rep(1:ncol(mat), each = nrow(mat)),
  is_gap        = as.vector(mat == "-")
)

# ordering sequences on y-axis per species
seq_order <- names(cds_seqs)[order(species)]
heatmap_long$sequence_name <- factor(heatmap_long$sequence_name, levels = seq_order)

ggplot(heatmap_long, aes(x = position, y = sequence_name, fill = is_gap)) +
  geom_tile() +
  scale_fill_manual(
    values = c("FALSE" = "grey20", "TRUE" = "tomato"),
    labels = c("FALSE" = "Nucleotide", "TRUE" = "Gap"),
    name = NULL
  ) +
  labs(
    title = "Gap distribution per sequence along the aligned CDS",
    x = "Alignment position",
    y = NULL
  ) +
  p14_theme +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = "top"
  )



# ==============================================================================
# 3. Session Info
# ==============================================================================

sessionInfo()
