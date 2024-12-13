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
chr_string       <- "NC_0475"
pedigree.FN      <- "04_impute/pedigree.csv"

# Load data
# Pedigree
pedigree.df <- read.table(file = pedigree.FN, header = FALSE, sep = " ")
colnames(pedigree.df) <- c("Animal", "Sire", "Dam")
pedigree.df$Sex <- NA
pedigree.df$Sex[grep(pattern = "F$", x = pedigree.df$Animal, perl = T)] <- "F"
pedigree.df$Sex[grep(pattern = "M$", x = pedigree.df$Animal, perl = T)] <- "M"
pedigree.df$Sex[is.na(pedigree.df$Sex)]  <- "M" # add M if unknown
write.table(x = pedigree.df, file = "04_impute/fimpute/pedigree.csv", sep = " ", quote = F, row.names = F, col.names = T)
# Only one pedigree file is needed, and can be reused for each chr

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
head(map.df)

# NOTE: here Chr does not correspond to the actual chromosome names of the species, it just keeps
#  in the same order as the NCBI identifiers for the chromosomes

# Make individual maps per chr
chunk <- NULL; map.FN <- NULL
for(i in 1:length(unique(map.df$Chr))){
  
  # Separate into chr
  chunk <- map.df[map.df$Chr==i,]
  chunk$Chip1 <- seq(1:nrow(chunk))
  
  # Create filename
  map.FN <- paste0("04_impute/fimpute/map_chr_", i, ".txt")
  
  # Write out
  write.table(x = chunk, file = map.FN, sep = " ", quote = F, row.names = F)
  
}


#### Genotype file ####
# Needs to be done per chr
chrs_to_include <- unique(gsub(pattern = " .*", replacement = "", x = ai2.df$mname))
chrs_to_include

# Per chr
geno_chunk <- NULL; geno_chunk_no_mname.df <- NULL; genotypes_no_mname.FN <- NULL
for(i in 1:length(chrs_to_include)){
  
  print(paste0("Working on chr ", i))
  
  # select just the section of interest
  geno_chunk <- ai2.df[grep(pattern = chrs_to_include[i], x = ai2.df$mname), ]
  
  # keep everything but mname
  geno_chunk_no_mname.df <- geno_chunk[, grep(pattern = "mname", x = colnames(geno_chunk), invert = T)]
  geno_chunk_no_mname.df <- t(geno_chunk_no_mname.df) # transpose
  geno_chunk_no_mname.df[1:5,1:5]
  
  # Make filename
  genotypes_no_mname.FN <- paste0("04_impute/fimpute/genotypes_no_indname_chr_", i, ".txt")
  
  # Write out the genotypes section
  fwrite(x = geno_chunk_no_mname.df, file = genotypes_no_mname.FN, sep = "", quote = F, row.names = F, col.names = F)
  
}

# The following two are constant, and can be reused per chr
# Inds section
indnames.df <- as.data.frame(rownames(geno_chunk_no_mname.df))
colnames(indnames.df) <- "ID"
head(indnames.df)
write.table(x = indnames.df, file = "04_impute/fimpute/indnames.txt", sep = " ", quote = F, row.names = F)

# Chip section
chip.df <- as.data.frame(rep(1, times = nrow(indnames.df)))
colnames(chip.df) <- "Chip"
head(chip.df)
write.table(x = chip.df, file = "04_impute/fimpute/chip.txt", sep = " ", quote = F, row.names = F)


# All component parts are prepared, now go to terminal and run 01_scripts/prep_fi3_indiv_chr.sh
