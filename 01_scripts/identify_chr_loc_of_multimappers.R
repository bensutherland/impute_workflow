# Identify the chr locations of multimappers from Sutherland et al. 2024 (File S4)
# B. Sutherland (2024-10-17)

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
bowtie_multimappers.FN <- "bowtie_multimappers.csv"
bwa_multimappers.FN    <- "bwa_multimappers.csv"
chr_pos_info.FN        <- "chr_pos_orig.mname.txt"

#### 01. Load data ####
# Read in multimappers
bowtie.df <- read.delim(file = bowtie_multimappers.FN, header = T, sep = ",")
head(bowtie.df)
nrow(bowtie.df)

bwa.df <- read.delim(file = bwa_multimappers.FN, header = T, sep = ",")
head(bwa.df)
nrow(bwa.df)

# Combine multimappers into a single vector
multimapper_mnames.vec <- c(bowtie.df$Marker.ID, bwa.df$Marker.ID)
length(multimapper_mnames.vec)

# Read in chr pos info
chr_pos_info.df <- read.delim(file = chr_pos_info.FN, header = F, sep = "\t")
head(chr_pos_info.df)
colnames(chr_pos_info.df) <- c("chr", "pos", "mname")

# Isolate the chr pos info to just those multimappers
caution_zones.df <- chr_pos_info.df[chr_pos_info.df$mname %in% multimapper_mnames.vec, ] 
nrow(caution_zones.df)

head(caution_zones.df)

write.table(x = caution_zones.df, file = "chr_locs_w_multimapper_loci.txt", sep = "\t", row.names = F)

