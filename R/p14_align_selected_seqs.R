# ==============================================================================
## ALIGNMENT OF SELECTED p14 CANDIDATES: METHOD COMPARISON ##
# ==============================================================================

# This script compares 5 alignment methods (DECIPHER, DECIPHER codon-aware, ClustalW,
# ClustalOmega, MUSCLE) on p14/ARF candidate CDS from several species, using
# length, gaps, conservation score and frame-consistency of indels.
# Plots each alignment against the human reference with exon 1b marked,
# compares nucleotide vs amino-acid (codon-level) conservation on
# the codon-aware alignment, and looks for conserved regions and
# whether they overlap exon 1b, as a simple check of support for homology.


library(Biostrings)
library(DECIPHER)
library(msa)
library(tidyverse)

# ==============================================================================
# 1. DATA LOADING
# ==============================================================================

fasta_dir  <- path.expand("~/Documents/PhD_work/CDKN2A_evo/raw_data/sequences/selected_p14_candidates")
output_dir <- path.expand("~/Documents/PhD_work/CDKN2A_evo/processed_data/alignments")

cds_seqs <- readDNAStringSet(file.path(fasta_dir, "p14_selected_species_combined_CDS.fasta"))

# Shorter names (the headers are the file names): "p14_candidate_mouse_CDS.fasta" -> "mouse"
names(cds_seqs) <- gsub("^p14_candidate_|_p14|_CDS|\\.fasta", "", names(cds_seqs))

# Human reference and exon 1b: aa 1-64 of human ARF are encoded by exon 1beta (192 nt)
ref_id    <- grep("human", names(cds_seqs), value = TRUE)
exon1b_nt <- 64 * 3


# test: remove anolis sequence
# cds_seqs <- cds_seqs[!grepl("anolis", names(cds_seqs), ignore.case = TRUE)]

# ==============================================================================
# 2. ALIGNMENTS
# ==============================================================================

aligned_decipher       <- AlignSeqs(cds_seqs, verbose = FALSE)
aligned_decipher_codon <- AlignTranslation(cds_seqs)   # translates, aligns the aa, back-translates
aligned_clustalw       <- msa(cds_seqs, method = "ClustalW")
aligned_clustalomega   <- msa(cds_seqs, method = "ClustalOmega")
aligned_muscle         <- msa(cds_seqs, method = "Muscle")

alignments <- list(
  DECIPHER       = aligned_decipher,
  DECIPHER_codon = aligned_decipher_codon,
  ClustalW       = DNAStringSet(aligned_clustalw),
  ClustalOmega   = DNAStringSet(aligned_clustalomega),
  MUSCLE         = DNAStringSet(aligned_muscle)
)

# writeXStringSet(aligned_decipher_codon, file.path(output_dir, "p14_alignment_decipher_codon.fasta"))
# writeXStringSet(aligned_decipher, file.path(output_dir, "p14_alignment_decipher.fasta"))

# ==============================================================================
# 3. COMPARING ALIGNMENT METHODS
# ==============================================================================

## 3.1 Basic statistics 

alignment_lengths <- data.frame(
  Method           = names(alignments),
  Alignment_length = sapply(alignments, function(x) width(x)[1]),
  row.names        = NULL
)

count_gaps <- function(alignment) {
  sum(sapply(as.character(alignment), function(x) str_count(x, "-")))
}

gap_percentage <- function(alignment) {
  total_gaps  <- count_gaps(alignment)
  total_bases <- sum(width(alignment))
  100 * total_gaps / total_bases
}

gap_stats <- data.frame(
  Method         = names(alignments),
  Gap_count      = sapply(alignments, count_gaps),
  Gap_percentage = sapply(alignments, gap_percentage),
  row.names      = NULL
)

## 3.2 Conservation 

alignment_quality <- function(alignment) {
  seq_matrix <- as.matrix(alignment)
  n_seqs     <- nrow(seq_matrix)

  quality <- data.frame(
    Position     = 1:ncol(seq_matrix),
    Conservation = NA,
    Occupancy    = NA
  )

  for (i in 1:ncol(seq_matrix)) {
    column       <- seq_matrix[, i]
    nucleotides  <- column[column != "-"]
    quality$Occupancy[i] <- length(nucleotides) / n_seqs
    if (length(nucleotides) > 0) {
      quality$Conservation[i] <- max(table(nucleotides)) / length(nucleotides)
    }
  }

  quality$Conservation_score <- quality$Conservation * quality$Occupancy * 100
  return(quality)
}

quality_list <- lapply(alignments, alignment_quality)

quality_summary <- data.frame(
  Method             = names(alignments),
  Mean_conservation  = sapply(quality_list, function(q) mean(q$Conservation, na.rm = TRUE)),
  Mean_occupancy     = sapply(quality_list, function(q) mean(q$Occupancy, na.rm = TRUE)),
  Mean_cons_score    = sapply(quality_list, function(q) mean(q$Conservation_score, na.rm = TRUE)),
  row.names          = NULL
)


## 3.3 Frame consistency of the indels 
# Only internal gaps count (gaps at the ends are not indels).
# With DECIPHER_codon this is ~100% by construction.

check_codon_frame <- function(alignment) {
  gap_lengths <- unlist(lapply(as.character(alignment), function(sequence) {
    runs        <- rle(strsplit(sequence, "")[[1]] == "-")
    n           <- length(runs$lengths)
    is_internal <- runs$values & seq_len(n) > 1 & seq_len(n) < n
    runs$lengths[is_internal]
  }))

  data.frame(
    Total_indels       = length(gap_lengths),
    Frame_consistent   = sum(gap_lengths %% 3 == 0),
    Frame_inconsistent = sum(gap_lengths %% 3 != 0)
  )
}

codon_frame_stats <- bind_rows(lapply(alignments, check_codon_frame), .id = "Method")
codon_frame_stats$Percent_consistent <- 100 * codon_frame_stats$Frame_consistent / codon_frame_stats$Total_indels

## 3.4 Summary table

alignment_summary <- alignment_lengths %>%
  left_join(gap_stats, by = "Method") %>%
  left_join(quality_summary, by = "Method") %>%
  left_join(codon_frame_stats, by = "Method")

#write.csv2(alignment_summary, file.path(output_dir, "alignment_summary.csv"), row.names = FALSE)

## 3.5 Plots

alignment_plot <- alignment_summary %>%
  select(Method, Alignment_length, Gap_count, Gap_percentage,
         Mean_cons_score, Total_indels, Percent_consistent) %>%
  pivot_longer(cols = -Method, names_to = "Metric", values_to = "Value") %>%
  mutate(Method = factor(Method, levels = names(alignments)))

p_alignment <- ggplot(alignment_plot, aes(x = Method, y = Value, fill = Method)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = round(Value, 1)), vjust = -0.35, size = 3) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 3) +
  scale_fill_viridis_d(option = "D") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Comparison of alignment methods", x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    strip.text         = element_text(face = "bold"),
    plot.title         = element_text(face = "bold", hjust = 0.5)
  )

#ggsave(filename = file.path(output_dir, "alignment_comparison.png"), plot = p_alignment, width = 12, height = 8, dpi = 300)

# Positional conservation score
positional_quality <- bind_rows(quality_list, .id = "Method") %>%
  mutate(Method = factor(Method, levels = names(alignments)))

p_position_alignment <- ggplot(
  positional_quality, aes(x = Position, y = Conservation_score, colour = Method)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ Method, ncol = 2) +
  scale_colour_viridis_d(option = "D") +
  scale_y_continuous(limits = c(0, 100)) +
  labs(title = "Positional conservation score", x = "Alignment position", y = "Conservation score") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text      = element_text(face = "bold"),
    plot.title      = element_text(face = "bold", hjust = 0.5)
  )

#ggsave(filename = file.path(output_dir, "positional_alignment_comparison.png"), plot = p_position_alignment, width = 12, height = 8, dpi = 300)

# ==============================================================================
# 4. ALIGNMENTS WITH EXON 1b MARKED
# ==============================================================================
# One plot per method, one tile per base, coloured by comparison with human.
# The black box marks the alignment columns covered by the first 192 nt (exon 1b)
# of the human sequence. 

n_seqs <- length(cds_seqs)

plot_alignment <- function(method) {
  m   <- as.matrix(alignments[[method]])
  ref <- m[ref_id, ]
  
  df <- expand.grid(Sequence = rownames(m), Position = 1:ncol(m), stringsAsFactors = FALSE)
  df$Base   <- m[cbind(match(df$Sequence, rownames(m)), df$Position)]
  df$Ref    <- ref[df$Position]
  
  df$Status <- "Different from human"
  df$Status[df$Base == df$Ref]     <- "Identical to human"
  df$Status[df$Sequence == ref_id] <- "Human (reference)"
  df$Status[df$Base == "-"]        <- "Gap"
  
  df$Status   <- factor(df$Status, levels = c("Human (reference)", "Identical to human",
                                              "Different from human", "Gap"))
  df$Sequence <- factor(df$Sequence, levels = rev(c(ref_id, setdiff(names(cds_seqs), ref_id))))
  
  # Columns of this alignment covered by exon 1b of human
  ref_cols <- which(ref != "-")
  xmin <- ref_cols[1] - 0.5
  xmax <- ref_cols[exon1b_nt] + 0.5
  
  ggplot(df, aes(x = Position, y = Sequence, fill = Status)) +
    geom_raster() +
    annotate("rect", xmin = xmin, xmax = xmax, ymin = 0.5, ymax = n_seqs + 0.5,
             fill = NA, colour = "black", linewidth = 0.8) +
    annotate("text", x = xmin, y = n_seqs + 0.6, hjust = 0, vjust = 0,
             label = "exon 1b (human aa 1-64)", fontface = "bold", size = 4) +
    scale_fill_manual(values = c("Human (reference)"    = "darkblue",
                                 "Identical to human"   = "lightgreen",
                                 "Different from human" = "lightcoral",
                                 "Gap"                  = "lightgrey")) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = expansion(add = c(0.5, 1.5))) +
    coord_cartesian(clip = "off") +
    labs(title = paste("Alignment:", method), x = "Alignment position", y = NULL, fill = NULL) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid      = element_blank(),
      plot.title      = element_text(face = "bold", hjust = 0.5),
      legend.position = "bottom"
    )
}

for (method in names(alignments)) {
  p <- plot_alignment(method)
  print(p)
 # ggsave(file.path(output_dir, paste0("alignment_exon1b_", method, ".png")),
 #       p, width = 14, height = 4, dpi = 300)
}

# ==============================================================================
# 5. CODON-LEVEL CONSERVATION (DECIPHER codon-aware alignment)
# ==============================================================================
# A codon is compared through the amino acid it encodes: a different nucleotide
# that gives the same amino acid still counts as conserved.
# Only works on the codon-aware alignment (its gaps are in triplets).

## 5.1 Codons and amino acids

nt_matrix <- as.matrix(aligned_decipher_codon)

# One column per codon (e.g. "ATG", "AAA", "---")
codon_matrix <- sapply(1:(ncol(nt_matrix) / 3), function(j) {
  apply(nt_matrix[, (3 * j - 2):(3 * j)], 1, paste, collapse = "")
})

# Amino acid encoded by each codon ("-" for gap codons)
aa_matrix <- matrix(GENETIC_CODE[codon_matrix], nrow = nrow(codon_matrix), dimnames = dimnames(codon_matrix))
aa_matrix[codon_matrix == "---"] <- "-"
aa_matrix[is.na(aa_matrix)]      <- "-"   # codons with ambiguous bases are treated as gaps

## 5.2 Conservation score (same function as for nucleotides, now on amino acids) 

aa_quality <- alignment_quality(aa_matrix)

data.frame(
  Level           = c("Nucleotide", "Codon (amino acid)"),
  Mean_cons_score = c(mean(quality_list$DECIPHER_codon$Conservation_score, na.rm = TRUE),
                      mean(aa_quality$Conservation_score, na.rm = TRUE))
)

# Last codon column covered by exon 1b of human
ref_cols   <- which(aa_matrix[ref_id, ] != "-")
exon1b_end <- ref_cols[exon1b_nt / 3]

codon_conservation_plot <- ggplot(aa_quality, aes(x = Position, y = Conservation_score)) +
  geom_line(linewidth = 0.7, colour = "steelblue") +
  geom_vline(xintercept = exon1b_end + 0.5, linetype = "dashed") +
  scale_y_continuous(limits = c(0, 100)) +
  labs(title = "Codon-level conservation score (DECIPHER codon-aware)",
       subtitle = "Dashed line = end of exon 1b",
       x = "Codon position in alignment", y = "Conservation score") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

#ggsave(filename = file.path(output_dir, "codon_conservation_score.png"), plot = codon_conservation_plot, width = 12, height = 8, dpi = 300)

## 5.3 Alignment plot by codon 

df <- expand.grid(Sequence = rownames(aa_matrix), Codon = 1:ncol(aa_matrix), stringsAsFactors = FALSE)
idx <- cbind(match(df$Sequence, rownames(aa_matrix)), df$Codon)

df$Codon_seq <- codon_matrix[idx]
df$Aa        <- aa_matrix[idx]
df$Ref_codon <- codon_matrix[ref_id, ][df$Codon]
df$Ref_aa    <- aa_matrix[ref_id, ][df$Codon]

df$Status <- "Different amino acid"
df$Status[df$Aa == df$Ref_aa]           <- "Same amino acid, different codon"
df$Status[df$Codon_seq == df$Ref_codon] <- "Identical codon"
df$Status[df$Sequence == ref_id]        <- "Human (reference)"
df$Status[df$Aa == "-"]                 <- "Gap"

df$Status   <- factor(df$Status, levels = c("Human (reference)", "Identical codon",
                                            "Same amino acid, different codon",
                                            "Different amino acid", "Gap"))
df$Sequence <- factor(df$Sequence, levels = rev(c(ref_id, setdiff(names(cds_seqs), ref_id))))

ggplot(df, aes(x = Codon, y = Sequence, fill = Status)) +
  geom_raster() +
  annotate("rect", xmin = ref_cols[1] - 0.5, xmax = exon1b_end + 0.5, ymin = 0.5, ymax = n_seqs + 0.5,
           fill = NA, colour = "black", linewidth = 0.8) +
  annotate("text", x = ref_cols[1] - 0.5, y = n_seqs + 0.6, hjust = 0, vjust = 0,
           label = "exon 1b (human aa 1-64)", fontface = "bold", size = 4) +
  scale_fill_manual(values = c("Human (reference)"                = "darkblue",
                               "Identical codon"                  = "lightgreen",
                               "Same amino acid, different codon" = "gold",
                               "Different amino acid"             = "lightcoral",
                               "Gap"                              = "lightgrey")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = expansion(add = c(0.5, 1.5))) +
  coord_cartesian(clip = "off") +
  labs(title = "DECIPHER codon-aware alignment compared with human, by codon",
       x = "Codon position in alignment", y = NULL, fill = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )

# ggsave(file.path(output_dir, "alignment_codons_decipher.png"), width = 14, height = 4, dpi = 300)

# ==============================================================================
# 6. CONSERVED REGIONS AND SUPPORT FOR HOMOLOGY
# ==============================================================================
# Uses the conservation score from section 3.2 (DECIPHER, nucleotides) and from
# section 5.2 (DECIPHER_codon, amino acids). Finds contiguous alignment positions
# with a score >= threshold, and checks whether they fall inside exon 1b.

threshold <- 70 # user-defined threshold 

find_conserved_regions <- function(quality, threshold) {
  above  <- replace(quality$Conservation_score >= threshold, is.na(quality$Conservation_score), FALSE)
  runs   <- rle(above)
  ends   <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1
  keep   <- runs$values   # only the TRUE runs, i.e. the conserved regions
  
  data.frame(
    Start      = starts[keep],
    End        = ends[keep],
    Length     = runs$lengths[keep],
    Mean_score = sapply(seq_along(starts[keep]), function(i) {
      mean(quality$Conservation_score[starts[keep][i]:ends[keep][i]])
    })
  )
}

nt_regions <- find_conserved_regions(quality_list$DECIPHER, threshold)
aa_regions <- find_conserved_regions(aa_quality, threshold)


## 6.1 Do the conserved regions fall inside exon 1b? 
# Exon 1b = alignment columns spanned by the first 192 nt (64 aa) of human.

nt_ref_cols <- which(as.matrix(aligned_decipher)[ref_id, ] != "-")
nt_exon1b   <- c(nt_ref_cols[1], nt_ref_cols[exon1b_nt])   # first, last column
aa_exon1b   <- c(ref_cols[1], exon1b_end)

overlaps_exon1b <- function(regions, exon1b) {
  regions$In_exon1b <- regions$Start <= exon1b[2] & regions$End >= exon1b[1]
  regions
}

nt_regions <- overlaps_exon1b(nt_regions, nt_exon1b)
aa_regions <- overlaps_exon1b(aa_regions, aa_exon1b)


## 6.2 Overall support for homology 
# How much of the alignment is well conserved, and how much of that sits in exon 1b.

homology_support <- data.frame(
  Level                     = c("Nucleotide (DECIPHER)", "Amino acid (DECIPHER codon-aware)"),
  Alignment_length          = c(nrow(quality_list$DECIPHER), nrow(aa_quality)),
  Mean_conservation_score   = c(mean(quality_list$DECIPHER$Conservation_score, na.rm = TRUE),
                                mean(aa_quality$Conservation_score, na.rm = TRUE)),
  Pct_positions_conserved   = c(100 * mean(quality_list$DECIPHER$Conservation_score >= threshold, na.rm = TRUE),
                                100 * mean(aa_quality$Conservation_score >= threshold, na.rm = TRUE)),
  N_conserved_regions       = c(nrow(nt_regions), nrow(aa_regions)),
  Conserved_positions_total = c(sum(nt_regions$Length), sum(aa_regions$Length)),
  Pct_conserved_in_exon1b   = c(100 * sum(nt_regions$Length[nt_regions$In_exon1b]) / sum(nt_regions$Length),
                                100 * sum(aa_regions$Length[aa_regions$In_exon1b]) / sum(aa_regions$Length))
)

#write.csv2(homology_support, file.path(output_dir, "homology_support_arf_seqs.csv"), row.names = FALSE)

## 6.3 Plot: conservation score with conserved regions and exon 1b highlighted 

plot_conservation <- function(quality, regions, exon1b, title) {
  ggplot(quality, aes(x = Position, y = Conservation_score)) +
    annotate("rect", xmin = exon1b[1] - 0.5, xmax = exon1b[2] + 0.5, ymin = -Inf, ymax = Inf,
             fill = "grey85") +
    geom_hline(yintercept = threshold, linetype = "dashed", colour = "grey40") +
    geom_line(colour = "steelblue") +
    geom_rect(data = regions, inherit.aes = FALSE, aes(xmin = Start - 0.5, xmax = End + 0.5),
              ymin = -3, ymax = 0, fill = "darkred") +
    annotate("text", x = mean(exon1b), y = 103, label = "human exon 1b", fontface = "bold", size = 4) +
    scale_y_continuous(limits = c(-3, 105)) +
    labs(title = title,
         subtitle = paste0("Dashed line = ", threshold, "% threshold; red bars = conserved regions"),
         x = "Alignment position", y = "Conservation score") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
}

nt_conservation_plot <- plot_conservation(
  quality_list$DECIPHER,
  nt_regions,
  nt_exon1b,
  "Nucleotide-level conservation (DECIPHER)"
)

aa_conservation_plot <- plot_conservation(
  aa_quality,
  aa_regions,
  aa_exon1b,
  "Amino-acid-level conservation (DECIPHER codon-aware)"
)

print(nt_conservation_plot)
print(aa_conservation_plot)

#ggsave(filename = file.path(output_dir, "nucleotide_conservation_DECIPHER.png"), plot = nt_conservation_plot, width = 12, height = 8, dpi = 300)
#ggsave(filename = file.path(output_dir, "amino_acid_conservation_DECIPHER.png"), plot = aa_conservation_plot, width = 12, height = 8, dpi = 300)


# ==============================================================================
# 7. SESSION INFO
# ==============================================================================

sessionInfo()
