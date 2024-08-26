# Read in ai2 file and convert to FImputev.3 file
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
ai2_input.FN     <- "../FImpute3_test/test_ai2_input_one_chr_subset.txt" 
ai2.df <- fread(file =  ai2_input.FN, sep = "\t")
dim(ai2.df)
ai2.df <- as.data.frame(ai2.df)
ai2.df[1:5,1:5]

# Map file
SNP <- gsub(pattern = " ", replacement = "__", x = ai2.df$mname)
Chr <- gsub(pattern = "\ .*", replacement = "", x = ai2.df$mname)
Pos <- gsub(pattern = ".*\ ", replacement = "", x = ai2.df$mname)

map.df <- cbind(SNP, Chr, Pos)
map.df <- as.data.frame(map.df)
map.df$Chip1 <- seq(1:nrow(map.df))
head(map.df)

write.table(x = map.df, file = "../FImpute3_test/map.txt", sep = " ", quote = F)

# Genotype file
ai2_no_mname.df <- ai2.df[, grep(pattern = "mname", x = colnames(ai2.df), invert = T)]
ai2_no_mname.df[1:5,1:5]
ai2_no_mname.df <- t(ai2_no_mname.df)
ai2_no_mname.df[1:5,1:5]


write.table(x = ai2_no_mname.df, file = "../FImpute3_test/genotypes_no_indname.txt", sep = "", quote = F
            , row.names = F, col.names = F
            )

indnames.df <- as.data.frame(rownames(ai2_no_mname.df))
colnames(indnames.df) <- "ID"
head(indnames.df)
write.table(x = indnames.df, file = "../FImpute3_test/indnames.txt", sep = " ", quote = F, row.names = F)

chip.df <- as.data.frame(rep(1, times = nrow(indnames.df)))
colnames(chip.df) <- "Chip"
head(chip.df)
write.table(x = chip.df, file = "../FImpute3_test/chip.txt", sep = " ", quote = F, row.names = F)
