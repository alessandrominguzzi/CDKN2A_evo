# ==============================================================================
## GFF PREP FOR MCScanX ##
# ==============================================================================

# This script aims to process the GFF3 files extracted from each species of 
# interest, in order to create a unique file suitable as an input for the 
# MCScanX software to perform synteny analysis. 

library(dplyr) 
library(tidyverse)

# ==============================================================================
# 1. DATA LOADING
# ==============================================================================

gff_dir <- file.path("~/Documents/PhD_work/CDKN2A_evo/raw_data/synteny_data/species_data_chr_specific/gff_files")
output_dir <- file.path("~/Documents/PhD_work/CDKN2A_evo/processed_data/synteny_mcscanx")

# checking files 
gff_files <- list.files(
  gff_dir,
  pattern = "\\.gff$",
  full.names = TRUE
)

# ==============================================================================
# 2. PROCESSING
# ==============================================================================

# Function to extract gene coordinates from GFFs
convert_gff <- function(file) {
  gff <- read.delim(
    file,
    header = FALSE,
    sep = "\t",
    comment.char = "#",
    quote = "",
    stringsAsFactors = FALSE
  )
  
  genes <- gff[gff$V3 == "gene", ] #keep rows with genes
  gene_ids <- sub(".*ID=gene-?([^;]+).*", "\\1", genes$V9) #extract ID
  
  data.frame(                #crate df format readable by MCScanX
    chromosome = genes$V1,
    gene_id = gene_ids,
    start = genes$V4,
    end = genes$V5
  )
}

# Unify all genes from all GFF
all_genes <- do.call(rbind, lapply(gff_files, convert_gff))

# Changing chr names 
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

all_genes <- all_genes %>%
  mutate(chromosome = unname(chromosome_map[chromosome]))

# Order by chr name
all_genes <- all_genes[order(all_genes$chromosome, all_genes$start), ]

#write.table(all_genes,file = file.path(output_dir, "all_species.gff"),sep = "\t", quote = FALSE,row.names = FALSE, col.names = FALSE)
