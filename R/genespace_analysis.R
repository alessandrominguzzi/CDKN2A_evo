# ==============================================================================
## GENESPACE ANALYSIS ##
# ==============================================================================
 
# This script is designed to perform an analysis of the genomic 
# coordinates of the identified candidate p16 sequences across different species. 
# The methodological approach is based on the R package GENESPACE.


library(Biostrings)
library(GENESPACE)
library(tidyverse)

# ==============================================================================
# 1. DATA LOADING
# ==============================================================================

## 1.1 Project directories

genespace_dir <- path.expand(
  "~/Documents/PhD_work/CDKN2A_evo/raw_data/genespace_data"
)

output_dir <- path.expand(
  "~/Documents/PhD_work/CDKN2A_evo/processed_data/synteny"
)

species_dirs <- basename(
  list.dirs(genespace_dir, recursive = FALSE)
)


## 1.2 Software directories

# Conda environment containing the GENESPACE dependencies
conda_dir <- path.expand(
  "~/miniforge3/envs/genespace"
)

# DIAMOND
diamond_path <- file.path(
  conda_dir, "bin", "diamond"
)

# OrthoFinder
orthofinder_path <- file.path(
  conda_dir, "bin", "orthofinder"
)

# MCScanX installation directory
mcscanx_path <- path.expand(
  "~/software/MCScanX"
)

conda_bin <- path.expand(
  "~/miniforge3/envs/genespace/bin"
)

Sys.setenv(
  PATH = paste(
    conda_bin,
    Sys.getenv("PATH"),
    sep = .Platform$path.sep
  )
)

orthofinder_path <- file.path(
  conda_bin, "orthofinder"
)

Sys.which("orthofinder")

system2(
  orthofinder_path,
  args = "-h",
  stdout = TRUE,
  stderr = TRUE
)


# ==============================================================================
# 2. DATA PROCESSING
# ==============================================================================

ParsedPath <- parse_annotations(
  rawGenomeRepo = genespace_dir,
  genomeDirs = species_dirs,
  genomeIDs = species_dirs,
  gffString = "gff",
  faString = "faa",
  presets = "ncbi",
  genespaceWd = output_dir
)

bed_dir <- file.path(
  output_dir, "bed"
)

peptide_dir <- file.path(
  output_dir, "peptide"
)


# ==============================================================================
# 3. GENESPACE
# ==============================================================================

## Initialize GENESPACE

gpar <- init_genespace(
  wd = output_dir,
  genomeIDs = species_dirs,
  path2diamond = diamond_path,
  path2orthofinder = orthofinder_path,
  path2mcscanx = mcscanx_path,
  nCores = 1
)

## Run GENESPACE

out <- run_genespace(gsParam = gpar)

# step 8 pan-genomes is not working
# step 4 (flagging synteny) required modified function for no synteny found for 
# some pairs of genomes 

# Orthofinder step OK





