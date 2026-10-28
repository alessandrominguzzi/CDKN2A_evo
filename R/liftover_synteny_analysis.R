# ==============================================================================
## LiftOver ANALYSIS - EVALUATING SYNTENY ##
# ==============================================================================

# This script investigates human-to-species-of-interest LiftOver chain files 
# obtained from the UCSC LiftOver tool. The aim is to obtain target
# sequences corresponding to p16/ARF sequences and investigate their synteny.
# The resulting coordinates and synteny results will be compared with 
# BLAST outputs and correspondent synteny analyses to evaluate consistency 
# and/or differences.

library(rtracklayer)
library(tidyverse)

# ==============================================================================
# 1. DATA LOADING
# ==============================================================================

chain_dir <- file.path("~/Documents/PhD_work/CDKN2A_evo/raw_data/synteny_data/chains")
input_bed <- file.path("~/Documents/PhD_work/CDKN2A_evo/raw_data/synteny_data/input.bed")
output_dir <- file.path("~/Documents/PhD_work/CDKN2A_evo/processed_data/synteny_liftover")

# ==============================================================================
# 2. LiftOver ANALYSIS 
# ==============================================================================

# defining regions of interest

regions <- import(input_bed) #converts 0-based to 1-based 
regions$name <- c("CDKN2A", "p16_locus", "ARF_locus", "exon1alpha", "exon1beta", "exon2", "exon3")

chain_files <- list.files(chain_dir, pattern = "over.chain$", full.names = TRUE)

# executing lifotver per chain file 

results <- data.frame()
for (f in chain_files) {
  chain  <- import.chain(f)
  lifted <- liftOver(regions, chain)
  for (i in seq_along(regions)) {
    pieces <- lifted[[i]]
    if (length(pieces) == 0) {
      row <- data.frame(chain = basename(f), region = regions$name[i],
                        chr = NA, start = NA, end = NA, strand = NA,
                        mapped_bp = 0, fraction = 0)
    } else {
      row <- data.frame(chain = basename(f), region = regions$name[i],
                        chr = as.character(seqnames(pieces)),
                        start = start(pieces), end = end(pieces),
                        strand = as.character(strand(pieces)),
                        mapped_bp = sum(width(pieces)),
                        fraction = sum(width(pieces)) / width(regions[i]))
    }
    results <- rbind(results, row)
  }
}

# write.csv(results, file.path(output_dir, "liftover_results.csv"), row.names = FALSE)

# processing results 

results$piece_bp <- results$end - results$start + 1   # bp of each lifted piece

summary_results <- results %>%
  filter(!is.na(chr)) %>%
  group_by(chain, region, chr) %>%
  summarise(start = min(start), end = max(end),  # unify fragments per region          
            mapped_bp = sum(piece_bp), n_pieces = n(), .groups = "drop") %>% # sum mapped bp and fragments per region
  group_by(chain, region) %>%
  slice_max(mapped_bp, n = 1, with_ties = FALSE) %>%   # keep the chr with most mapped bp
  ungroup()

# one row per chain, regions side by side
summary_wide <- summary_results %>%
  pivot_wider(id_cols = chain, names_from = region,
              values_from = c(chr, start, end))

# write.csv(summary_results, file.path(output_dir, "liftover_summary.csv"), row.names = FALSE)
# write.csv(summary_wide, file.path(output_dir, "liftover_summary_wide.csv"), row.names = FALSE)

# isolating CDKN2A, p16, ARF best results 

CDKN2A_res <- summary_wide %>% 
  select(chain, matches("CDKN2A"))

ARF_res <- summary_wide %>% 
  select(chain, matches("ARF"))

p16_res <- summary_wide %>% 
  select(chain, matches("p16"))
