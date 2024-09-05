# Prepare imputed ai2-formatted file for GEMMA analysis
# B. Sutherland
# initialized 2024-07-26

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

# Set user variables
#imputed.FN     <- "04_impute/fimpute/fi3_loci_by_inds_all_imputed_chr.txt" # imputed file in ai2 format
imputed.FN     <- "05_compare/all_chr_combined.txt" # imputed file in ai2 format
phenotypes.FN <- "00_archive/G0923-21-VIUN_SampleInventory_V2_recd_2024-08-16.txt"
offspring_string  <- "ASY2"

pheno_of_interest <- "survival_state" # DPE or survival_state


#### 01. Load data ####
## Read in imputed data
imputed.df <- fread(file = imputed.FN, sep = "\t")
dim(imputed.df)
imputed.df <- as.data.frame(imputed.df) # convert to df
imputed.df[1:5,1:5]

# Format
colnames(imputed.df)[1] <- "mname"
imputed.df$mname <- gsub(pattern = " ", replacement = "__", x = imputed.df$mname) # will correct mname if needed

# Drop all samples but offspring (assumes first column is mname)
imputed.df <- imputed.df[ , c(1, grep(pattern = offspring_string, x = colnames(imputed.df)))]
dim(imputed.df)
imputed.df[1:5,1:5]

## Read in phenotypes
phenos.df <- read.delim(file = phenotypes.FN)
phenos.df$indiv <- gsub(pattern = "_", replacement = "-", x = phenos.df$indiv) # make matching to imputed colnames
phenos.df <- as.data.frame(phenos.df)
head(phenos.df)
nrow(phenos.df) # 240 inds


#### 02. Prepare GEMMA input files ####
##### a. genotype file #####
# Add dummy cols required by GEMMA
imputed.df$x <- "X"
imputed.df$y <- "Y"
colnames(imputed.df)

# Sort by colname, then put mname first
imputed.df <- imputed.df %>% 
  select("mname", "x", "y", everything())
imputed.df[1:5,1:5]

# Ensure all indivs have phenos before proceeding
missing_pheno_inds <- phenos.df$indiv[is.na(phenos.df[, pheno_of_interest])]

# Drop any indivs from the imputed file if they have missing data
imputed.df <- imputed.df[, !(colnames(imputed.df) %in% missing_pheno_inds)]
dim(imputed.df)

## Write output imputed genos
fwrite(x = imputed.df, file = "07_GWAS/gwas_geno.txt", sep = " ", col.names = F)

##### b. phenotype file #####
## Prepare pheno file
indivs.df <- colnames(imputed.df)[grep(pattern = "mname|x|y", x = colnames(imputed.df), invert = T)]
indivs.df <- as.data.frame(indivs.df)
colnames(indivs.df) <- "indiv"
head(indivs.df)

# Are all inds in the impute file present in the NA-removed pheno file? 
length(intersect(x = indivs.df$indiv, y = phenos.df$indiv))

# Add phenotypes to the indiv file (in the order they are in the imputed file)
ordered_phenos.df <- merge(x = indivs.df, y = phenos.df, by = "indiv", all.x = T, sort = F)
ordered_phenos.df

# Inspect briefly to ensure order is retained
head(ordered_phenos.df)
head(colnames(imputed.df), n = 10)
tail(ordered_phenos.df)
tail(colnames(imputed.df), n = 6)

# Update survivor to day 17
if(pheno_of_interest=="DPE"){
  
  print("Update survivors' DPE to a day after the trial")
  ordered_phenos.df[ordered_phenos.df$survival_state=="S", "DPE"] <- "17"
  var_status <- as.numeric(ordered_phenos.df$DPE)
  
}else if(pheno_of_interest=="survival_state"){
  
  ordered_phenos.df$survival_state <- gsub(pattern = "S", replacement = "1", x = ordered_phenos.df$survival_state)
  ordered_phenos.df$survival_state <- gsub(pattern = "M", replacement = "0", x = ordered_phenos.df$survival_state)
  var_status <- as.numeric(ordered_phenos.df$survival_state)
  
}

gwaspheno <-  var_status
write.table(x = gwaspheno, file = "07_GWAS/gwas_pheno.txt", row.names = F, col.names = F)


##### c. covariate file #####
# Obtain vector of family identities
inds_present <- ordered_phenos.df$indiv
inds_present[grep(pattern = "-114-", x = inds_present)] <- "F114"
inds_present[grep(pattern = "-115-", x = inds_present)] <- "F115"
inds_present[grep(pattern = "-116-", x = inds_present)] <- "F116"
inds_present[grep(pattern = "-117-", x = inds_present)] <- "F117"
table(inds_present)

# Retain as covariate
gwascovar = model.matrix(~as.factor(inds_present))
write.table(x = gwascovar, file = "07_GWAS/gwas_covar.txt", row.names = F, col.names = F)

# Now run gemma
