# Convert AI2 output back to VCF format (final step)
# note: requires that ./01_scripts/ai2_to_vcf_format_step_1.sh has been run, and that 
#  a VCF file with only the imputed loci has been generated
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
## AI2 settings 
# pre_impute_VCF.FN  <- "04_impute/all_inds_wgrs_and_panel_biallele_only_ai2_imputed_regions.vcf" # this is the pre-impute VCF, subset to only incl the loci retained after the ai2 imputation
# post_impute_ai2.FN <- "05_compare/all_chr_combined_converted.txt"

## FI3 settings
pre_impute_VCF.FN  <- "04_impute/all_inds_wgrs_and_panel_biallele_only_fi3_imputed_regions.vcf" # this is the pre-impute VCF, subset to only incl the loci retained after the fi3 imputation
post_impute_ai2.FN <- "04_impute/fimpute/fi3_loci_by_inds_all_imputed_chr_converted.txt"


#### 01. Load data ####
## AI2-format imputed data ##
# Load AI2 output file
ai2_output.df <- fread(file = post_impute_ai2.FN, sep = "\t")
dim(ai2_output.df)
ai2_output.df$mname <- gsub(pattern = " ", replacement = "__", x = ai2_output.df$mname) # to match the VCF names
ai2_output.df[1:5,1:5]
gc()

# Load VCF file
pre_vcf <- read.vcfR(file = pre_impute_VCF.FN)
gc()

# View GT section
dim(pre_vcf@gt)
pre_vcf@gt[1:5,1:5]

# NOTE: the AI2 and the pre-imputed VCF may NOT have the same sample order, as it may have been 
#  sorted by alphaimpute2 (2024-12-13, bjgs)

# View Info section
dim(pre_vcf@fix)
head(pre_vcf@fix)

# Extract marker order in pre-imputed VCF
vcf_genos.df <- paste0(pre_vcf@fix[, "CHROM"], "__", pre_vcf@fix[, "POS"])
vcf_genos.df <- as.data.frame(vcf_genos.df)
colnames(vcf_genos.df) <- "mname"
head(vcf_genos.df)
dim(vcf_genos.df)

# Merge GT with ai2 (must not sort, to put in order of input VCF)
replacement_genos.df <- merge(x = vcf_genos.df, y = ai2_output.df, by = "mname", all.x = T, sort = F)
dim(replacement_genos.df)
replacement_genos.df[1:5,1:5]

# Note: The matrix should be in the correct order for the markers/marker names in the VCF

# Make new 'Format' column (that says GT only)
# Replace the mname column with the FORMAT
colnames(replacement_genos.df) <- gsub(pattern = "mname", replacement = "FORMAT", x = colnames(replacement_genos.df))
replacement_genos.df$FORMAT <- rep("GT", times = nrow(replacement_genos.df))
replacement_genos.df[1:5,1:5]

# Put the columns in the same order as expected from the pre-impute VCF (ai2 shifted them)
replacement_genos.df <- replacement_genos.df[, colnames(pre_vcf@gt)]
replacement_genos.df[1:5,1:10]

# Replace the GT section of the VCF
dim(pre_vcf@gt)
dim(replacement_genos.df)
replacement_genos.mat <- as.matrix(replacement_genos.df)
pre_vcf@gt <- replacement_genos.mat

pre_vcf@gt[1:10,1:10]

# Write VCF
output.FN <- gsub(pattern = "_imputed_regions.vcf", replacement = "_imputed.vcf.gz", x = pre_impute_VCF.FN)
write.vcf(x = pre_vcf, file = output.FN)

# It is now possible to compare the VCF file to other VCF files
