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
pre_impute_VCF.FN  <- "04_impute/all_inds_wgrs_and_panel_biallele_only_ai2_imputed_regions.vcf" # this is the pre-impute VCF, subset to only incl the loci retained after the ai2 imputation
post_impute_ai2.FN <- "05_compare/all_chr_combined_converted.txt"


#### 01. Load data ####
# Load VCF file
pre_vcf <- read.vcfR(file = pre_impute_VCF.FN)
#pre_vcf.bck <- pre_vcf
gc()

# Load AI2 output file
ai2_output.df <- fread(file = post_impute_ai2.FN, sep = "\t")
dim(ai2_output.df)
ai2_output.df$mname <- gsub(pattern = " ", replacement = "__", x = ai2_output.df$mname) # to match the VCF names
ai2_output.df[1:5,1:5]
gc()

# View GT section
dim(pre_vcf@gt)
pre_vcf@gt[1:5,1:5]

# View Info section
dim(pre_vcf@fix)
head(pre_vcf@fix)

# Extract GT
#pre_vcf_gt.df <- extract.gt(x = pre_vcf, element = "GT") ### THIS MAY NOT BE NECESSARY, JUST USE THE ORDER!!! 
vcf_genos.df <- paste0(pre_vcf@fix[, "CHROM"], "__", pre_vcf@fix[, "POS"])
vcf_genos.df <- as.data.frame(vcf_genos.df)
colnames(vcf_genos.df) <- "mname"
head(vcf_genos.df)

# Merge GT with ai2
replacement_genos.df <- merge(x = vcf_genos.df, y = ai2_output.df, by = "mname", all.x = T, sort = F)

# Make new Format column, just GT
dim(replacement_genos.df)
replacement_genos.df[1:5,1:5]

# Replace the mname column with the FORMAT
colnames(replacement_genos.df) <- gsub(pattern = "mname", replacement = "FORMAT", x = colnames(replacement_genos.df))
replacement_genos.df$FORMAT <- rep("GT", times = nrow(replacement_genos.df))

# Replace the GT section of the VCF
dim(pre_vcf@gt)
dim(replacement_genos.df)
replacement_genos.mat <- as.matrix(replacement_genos.df)
pre_vcf@gt <- replacement_genos.mat

pre_vcf@gt[1:10,1:10]

# Write VCF
write.vcf(x = pre_vcf, file = "05_compare/all_inds_wgrs_and_panel_biallele_imputed.vcf.gz")
