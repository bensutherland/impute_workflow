# Convert AI2 output back to VCF format
# B. Sutherland (2024-09-04)


### Front Matter ####
# Clean space
# rm(list=ls())

## Install and load packages
#install.packages("rstudioapi")
library(rstudioapi)
library(vcfR)
library(data.table)

# Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "/01_scripts", replacement = "", x = current.path)
setwd(current.path)
rm(current.path)

# Set user variables
pre_impute_VCF.FN  <- "04_impute/all_inds_wgrs_and_panel_biallele.vcf"
post_impute_ai2.FN <- "05_compare/all_chr_combined.txt"

# Load VCF file
pre_vcf <- read.vcfR(file = pre_impute_VCF.FN)

# Extract GT section
dim(pre_vcf@gt)
pre_vcf@gt[1:5,1:5]
dim(pre_vcf@fix)
pre_vcf@fix[1:5,1:5]

# Load AI2 output file
ai2_output.df <- fread(file = post_impute_ai2.FN, sep = "\t")
dim(ai2_output.df)
ai2_output.df[1:5,1:5]
# mname is space sep




