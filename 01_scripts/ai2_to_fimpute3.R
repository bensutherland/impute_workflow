# Read in ai2 file and convert to FImpute3 file
# B. Sutherland (2024-08-26)

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

# Set variables
ai2_input.FN     <- "04_impute/all_inds_wgrs_and_panel_biallele_ai2.txt" 
chr_string <- "NC_0475"
pedigree.FN <- "04_impute/pedigree.csv"

# Load data
# Pedigree
pedigree.df <- read.table(file = pedigree.FN, header = FALSE, sep = " ")
colnames(pedigree.df) <- c("Animal", "Sire", "Dam")
pedigree.df$Sex <- NA
pedigree.df$Sex[grep(pattern = "F$", x = pedigree.df$Animal, perl = T)] <- "F"
pedigree.df$Sex[grep(pattern = "M$", x = pedigree.df$Animal, perl = T)] <- "M"
pedigree.df$Sex[is.na(pedigree.df$Sex)]  <- "M" # add M if unknown
write.table(x = pedigree.df, file = "04_impute/fimpute/pedigree.csv", sep = " ", quote = F, row.names = F, col.names = T)

# AI2
ai2.df <- fread(file =  ai2_input.FN, sep = "\t")
dim(ai2.df)
ai2.df <- as.data.frame(ai2.df)
ai2.df[1:5,1:5]

# Only keep chromosomes
ai2.df <- ai2.df[grep(pattern = chr_string, x = ai2.df$mname), ]
dim(ai2.df)

#### Prepare the map file ####
SNP <- gsub(pattern = " ", replacement = "__", x = ai2.df$mname)

# Chr must be numeric... 
Chr <- gsub(pattern = "\ .*", replacement = "", x = ai2.df$mname)
Chr <- gsub(pattern = "NC_047559.1", replacement = 1, x = Chr)
Chr <- gsub(pattern = "NC_047560.1", replacement = 2, x = Chr)
Chr <- gsub(pattern = "NC_047561.1", replacement = 3, x = Chr)
Chr <- gsub(pattern = "NC_047562.1", replacement = 4, x = Chr)
Chr <- gsub(pattern = "NC_047563.1", replacement = 5, x = Chr)
Chr <- gsub(pattern = "NC_047564.1", replacement = 6, x = Chr)
Chr <- gsub(pattern = "NC_047565.1", replacement = 7, x = Chr)
Chr <- gsub(pattern = "NC_047566.1", replacement = 8, x = Chr)
Chr <- gsub(pattern = "NC_047567.1", replacement = 9, x = Chr)
Chr <- gsub(pattern = "NC_047568.1", replacement = 10, x = Chr)
unique(Chr)

Pos <- gsub(pattern = ".*\ ", replacement = "", x = ai2.df$mname)

map.df <- cbind(SNP, Chr, Pos)
map.df <- as.data.frame(map.df)
map.df$Chip1 <- seq(1:nrow(map.df))
head(map.df)
dim(map.df)

write.table(x = map.df, file = "04_impute/fimpute/map.txt", sep = " ", quote = F, row.names = F)


#### Genotype file ####
# Genotypes section
ai2_no_mname.df <- ai2.df[, grep(pattern = "mname", x = colnames(ai2.df), invert = T)] # keep everything but mname
ai2_no_mname.df[1:5,1:5]
ai2_no_mname.df <- t(ai2_no_mname.df) # transpose
ai2_no_mname.df[1:5,1:5]

fwrite(x = ai2_no_mname.df, file = "04_impute/fimpute/genotypes_no_indname.txt", sep = "", quote = F, row.names = F, col.names = F)

# Inds section
indnames.df <- as.data.frame(rownames(ai2_no_mname.df))
colnames(indnames.df) <- "ID"
head(indnames.df)
write.table(x = indnames.df, file = "04_impute/fimpute/indnames.txt", sep = " ", quote = F, row.names = F)

# Chip section
chip.df <- as.data.frame(rep(1, times = nrow(indnames.df)))
colnames(chip.df) <- "Chip"
head(chip.df)
write.table(x = chip.df, file = "04_impute/fimpute/chip.txt", sep = " ", quote = F, row.names = F)


# All component parts are prepared, now go to terminal and run 01_scripts/prep_fi3.sh
