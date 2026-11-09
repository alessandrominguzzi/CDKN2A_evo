# ==============================================================================
## FILES PREP FOR MCScanX ##
# ==============================================================================

# This script aims to process the GFF3 and FASTA (protein) files extracted from each species of 
# interest, in order to create unique files suitable as inputs for the 
# MCScanX software to perform synteny analysis. 

library(dplyr)
library(tidyverse)
library(Biostrings)

# ==============================================================================
# 1. DATA LOADING
# ==============================================================================

# Define directories
gff_dir <- file.path("~/Documents/PhD_work/CDKN2A_evo/raw_data/synteny_data/species_data_chr_specific/gff_files")
fasta_prot_dir <- file.path("~/Documents/PhD_work/CDKN2A_evo/raw_data/synteny_data/species_data_chr_specific/proteins_files")
output_dir <- file.path("~/Documents/PhD_work/CDKN2A_evo/processed_data/synteny_mcscanx")

# List files
gff_files <- list.files(gff_dir, pattern = "\\.gff$", full.names = TRUE)
faa_files <- list.files(fasta_prot_dir, pattern = "\\.faa$", full.names = TRUE)

# ==============================================================================
# 2. PROCESSING GFF 
# ==============================================================================

# Define chromosome names
chromosome_map <- c(
  "NC_000009.12" = "Hsa_chr9",
  "NC_000070.7"  = "Mmu_chr4",
  "NC_030724.1"  = "Xla_chr1L",
  "NC_041753.1"  = "Oan_chrX5",
  "NC_045790.1"  = "Bmu_chr6",
  "NC_052572.1"  = "Gga_chrZ",
  "NC_077232.2"  = "Mdom_chr6",
  "NC_085842.1"  = "Acar_chr2",
  "NC_133176.1"  = "Dre_chr1",
  "NC_134716.1"  = "Psin_chr6"
)

# Define species prefixes
species_map <- sub("_chr.*", "", chromosome_map)

# Convert GFF files
convert_gff <- function(file) {
  
  gff <- read.delim(file, header = FALSE, sep = "\t", comment.char = "#",
                    quote = "", stringsAsFactors = FALSE)
  
  genes <- gff[gff$V3 == "gene", ]
  
  gene_ids <- sub(".*ID=gene-?([^;]+).*", "\\1", genes$V9)
  
  chromosome <- unname(chromosome_map[genes$V1])
  species <- unname(species_map[genes$V1])
  
  gene_ids <- paste(species, gene_ids, sep = "_")
  
  data.frame(
    chromosome = chromosome,
    gene_id = gene_ids,
    start = genes$V4,
    end = genes$V5
  )
}

# Combine all GFF files
all_genes <- do.call(rbind, lapply(gff_files, convert_gff))

# Sort genes by chromosome and position
all_genes <- all_genes %>%
  arrange(chromosome, start)

# ==============================================================================
# 3. PROCESSING PROTEINS 
# ==============================================================================

# Define species prefixes for FAA files
# Replace the filenames below with the exact names returned by basename(faa_files)
faa_species_map <- c(
  "human_prot_ch9.faa" = "Hsa",
  "mouse_prot_ch4.faa" = "Mmu",
  "xenopus_prot_ch1L.faa" = "Xla",
  "platypus_prot_chX5.faa" = "Oan",
  "blue_whale_prot_ch6.faa" = "Bmu",
  "chicken_prot_chZ.faa" = "Gga",
  "opossum_prot_ch6.faa" = "Mdom",
  "anolis_prot_ch2.faa" = "Acar",
  "zebrafish_prot_ch1.faa" = "Dre",
  "pelodiscus_prot_ch6.faa" = "Psin"
)

# Convert FAA files
convert_faa <- function(file) {
  
  proteins <- readAAStringSet(file)
  
  gene_ids <- str_match(names(proteins), "\\[gene=([^]]+)\\]")[, 2]
  
  # Remove proteins without gene IDs
  valid <- !is.na(gene_ids) & gene_ids != ""
  proteins <- proteins[valid]
  gene_ids <- gene_ids[valid]
  
  # Select the longest protein for each gene
  protein_lengths <- width(proteins)
  
  selected <- tibble(
    index = seq_along(proteins),
    gene_id = gene_ids,
    length = protein_lengths
  ) %>%
    group_by(gene_id) %>%
    slice_max(order_by = length, n = 1, with_ties = FALSE) %>%
    pull(index)
  
  proteins <- proteins[selected]
  gene_ids <- gene_ids[selected]
  
  # Add species prefix to gene IDs
  filename <- basename(file)
  species <- unname(faa_species_map[filename])
  gene_ids <- paste(species, gene_ids, sep = "_")
  names(proteins) <- gene_ids
  proteins
}

# Combine all protein files
all_proteins <- do.call(c, lapply(faa_files, convert_faa))

# ==============================================================================
# 4. QUALITY CONTROLS AND EXPORT
# ==============================================================================

# Check results
length(all_genes) 
length(all_proteins) #14928

# Check duplicated gene IDs
sum(duplicated(all_genes$gene_id))  #no duplicates
sum(duplicated(names(all_proteins))) #no duplicates 

# Check common gene IDs between GFF and FAA files
gff_ids <- all_genes$gene_id
faa_ids <- names(all_proteins)

length(intersect(gff_ids, faa_ids)) #14191
length(setdiff(gff_ids, faa_ids))   #6293
length(setdiff(faa_ids, gff_ids))  #9

# Check chromosome names
unique(all_genes$chromosome)

# genes in the GFF without proteins are kept. Proteins without associated genomic
# coordinates (genes) are filtered out.
common_ids <- intersect(
  all_genes$gene_id,
  names(all_proteins)
)

proteins_filtered <- all_proteins[
  names(all_proteins) %in% common_ids
]


# write.table(all_genes, file = file.path(output_dir, "all_species.gff"), sep = "\t",quote = FALSE, row.names = FALSE, col.names = FALSE)
# writeXStringSet(proteins_filtered, filepath = file.path(output_dir, "all_species.faa"))


