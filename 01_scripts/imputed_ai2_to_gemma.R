# Prepare imputed ai2 file for gemma
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

# Set variables
offspring_imputed_ai2.FN     <- "13_impute_compare/all_chr_combined.txt" # imputed
phenotypes.FN <- "00_archive/qcat992_sample_mort_pheno_2024-06-17.txt"
pheno_of_interest <- "days_post_challenge"
  

#### 01. Load data ####
# Read in imputed data
imputed.df <- fread(file = offspring_imputed_ai2.FN, sep = "\t")
dim(imputed.df)
imputed.df <- as.data.frame(imputed.df) # convert to df
imputed.df[1:5,1:5]


#### 02. Prepare data for matching ####
head(colnames(imputed.df), n = 20)

# Remove parents
imputed.df <- imputed.df[, grep(pattern = "mname|ASY2", x = colnames(imputed.df))] # remove parents, keep mname and ASY2 inds
dim(imputed.df)

# Remove string '_ReAMP' from the end of sample names
colnames(imputed.df) <- gsub(pattern = "_ReAMP", replacement = "", x = colnames(imputed.df))
head(colnames(imputed.df))
table(duplicated(colnames(imputed.df))) # any duplicates?

head(colnames(imputed.df), n = 20)

#### 04. Clean up mname
imputed.df$mname <- gsub(pattern = " ", replacement = "__", x = imputed.df$mname)
imputed.df[1:5,1:5]


#### 05. Prepare for gemma ####
# Add dummy cols
imputed.df$x <- "X"
imputed.df$y <- "Y"
colnames(imputed.df)

# Sort by colname, then put mname first
imputed.df <- imputed.df %>% 
  select("mname", "x", "y", everything())
imputed.df[1:5,1:5]

# Need to make sure all indivs have phenos before proceeding
# Also need to convert NA to 9
#fwrite(x = imputed.df, file = "14_imputed_gwas/imputed_geno.txt", sep = " ", col.names = F)


#### 06. Prepare covariate and pheno files
pheno.df <- colnames(imputed.df)[4:ncol(imputed.df)]
pheno.df <- as.data.frame(pheno.df)

true_phenos.df <- read.delim(file = phenotypes.FN)
head(true_phenos.df)
nrow(true_phenos.df)

length(intersect(x = pheno.df$pheno.df, y = true_phenos.df$tube_label))
# So some are missing
drop.inds <- setdiff(x = pheno.df$pheno.df, y = true_phenos.df$tube_label)
drop.inds

imputed_w_pheno.df <- imputed.df[, !(colnames(imputed.df) %in% drop.inds)]
colnames(imputed_w_pheno.df)

pheno.df <- colnames(imputed_w_pheno.df)[4:ncol(imputed_w_pheno.df)]
pheno.df <- as.data.frame(pheno.df)

setdiff(x = pheno.df$pheno.df, y = true_phenos.df$tube_label) # OK

pheno.df <- merge(x = pheno.df, y = true_phenos.df, by.x = "pheno.df", by.y = "tube_label", all.x = T, sort = F)
pheno.df
head(cbind(colnames(imputed_w_pheno.df)[4:ncol(imputed_w_pheno.df)], pheno.df))
tail(cbind(colnames(imputed_w_pheno.df)[4:ncol(imputed_w_pheno.df)], pheno.df))

# OK, in order

# Update survivor to day 17
pheno.df[pheno.df$mort_surv=="S", "days_post_challenge"] <- "17"
head(pheno.df)

if(pheno_of_interest=="mort_surv"){
  
  pheno.df$mort_surv <- gsub(pattern = "S", replacement = "1", x = pheno.df$mort_surv)
  pheno.df$mort_surv <- gsub(pattern = "M", replacement = "0", x = pheno.df$mort_surv)
  var_status <- as.numeric(pheno.df$mort_surv)
  var_status
  
}else if(pheno_of_interest=="days_post_challenge"){
  
  var_status <- as.numeric(pheno.df$days_post_challenge)
  
}

gwaspheno <-  var_status
write.table(x = gwaspheno, file = "14_imputed_gwas/gwas_pheno.txt", row.names = F, col.names = F)

## Obtain vector of family identities
inds_present <- colnames(imputed_w_pheno.df)[4:ncol(imputed_w_pheno.df)]
inds_present[grep(pattern = "_114_", x = inds_present)] <- "F114"
inds_present[grep(pattern = "_115_", x = inds_present)] <- "F115"
inds_present[grep(pattern = "_116_", x = inds_present)] <- "F116"
inds_present[grep(pattern = "_117_", x = inds_present)] <- "F117"
table(inds_present)

# Retain as covariate
gwascovar = model.matrix(~as.factor(inds_present))
write.table(x = gwascovar, file = "14_imputed_gwas/gwas_covar.txt", row.names = F, col.names = F)

# Also need to convert NA to 9 (Not needed, since there are no missing values after imputation)

fwrite(x = imputed_w_pheno.df, file = "14_imputed_gwas/imputed_geno.txt", sep = " ", col.names = F)


# Now run gemma