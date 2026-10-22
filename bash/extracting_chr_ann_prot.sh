{\rtf1\ansi\ansicpg1252\cocoartf2822
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 #!/bin/bash\
\
# Chromosome accessions are stored in the metadata (Excel files).\
# The following commands are examples only: human chromosome 9\
# for annotation extraction (gff file) and protein extraction (faa file).\
# The same procedure must be performed for every species included in the study.\
\
# Extract annotations for the chromosome of interest \
awk '$1 == "#" || $1 == "NC_000009.12"' human_genomic.gff > human_chr9.gff\
\
# Extract proteins for the chromosome of interest \
seqkit grep -r -p "NC_000009.12" human.faa > human_ch9.faa}