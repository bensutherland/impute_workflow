# Read in AlphaImpute input and drop specific loci
# B. Sutherland (2024-08-01)

### Front Matter ####
# Clean space
# rm(list=ls())

## Install and load packages
#install.packages("rstudioapi")
#install.packages("data.table")
library("rstudioapi")
library("data.table")
library(tidyr)
library(dplyr)

# Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "/01_scripts", replacement = "", x = current.path)
setwd(current.path)
rm(current.path)

## Info
# sessionInfo()

# Set variables
ai2_input.FN    <- "12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic_no_MERR_ai2.txt"
bowtie_multimappers.FN <- "12_impute_impute/bowtie_multimappers.csv"
bwa_multimappers.FN <- "12_impute_impute/bwa_multimappers.csv"
hotspot_loci_in_bcf.FN <- "12_impute_impute/hotspot_loci_in_bcf.txt"

#### 01. Load data ####
# Read in ai2 input file (with all loci)
ai2_input.df <- fread(file = ai2_input.FN, sep = "\t")
dim(ai2_input.df)
ai2_input.df <- as.data.frame(ai2_input.df) # convert to df
ai2_input.df[1:5,1:5]

# Read in multimappers
bowtie.df <- read.delim(file = bowtie_multimappers.FN, header = T, sep = ",")
head(bowtie.df)

bwa.df <- read.delim(file = bwa_multimappers.FN, header = T, sep = ",")
head(bwa.df)

# Combine multimappers into a single vector
multimapper_mnames.vec <- c(bowtie.df$Marker.ID, bwa.df$Marker.ID)

# Read in BCF of hotspot loci
hotspot_loci.df <- read.delim(file = hotspot_loci_in_bcf.FN, header = F, sep = "\t")
dim(hotspot_loci.df)
hotspot_loci.df <- hotspot_loci.df[,1:3]
head(hotspot_loci.df)

# Subset to only target the multimappers in the bcf hotspot info
hotspot_loci_to_rem.df <- hotspot_loci.df[hotspot_loci.df$V3 %in% multimapper_mnames.vec, ]
head(hotspot_loci_to_rem.df)
table(hotspot_loci_to_rem.df$V3 %in% multimapper_mnames.vec)
length(multimapper_mnames.vec)

# Create matching vector to drop from the ai2 file
hotspot_loci_to_rem_chr_pos.vec <- paste0(hotspot_loci_to_rem.df$V1, " ", hotspot_loci_to_rem.df$V2)

head(ai2_input.df)[1:5]
ai2_input_clean.df <- ai2_input.df[!(ai2_input.df$mname %in% hotspot_loci_to_rem_chr_pos.vec), ]
dim(ai2_input.df)
dim(ai2_input_clean.df)

# Write output
fwrite(x = ai2_input_clean.df, file = "12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic_no_MERR_no_MM_ai2.txt", sep = "\t", quote = F, col.names = T)
