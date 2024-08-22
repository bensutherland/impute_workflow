# Convert the geno matrix to inputs needed for AlphaImpute2
#  note: requires that 01_scripts/bcf_to_ai2.sh has already been run
# B. Sutherland (2024-07-19)

### Front Matter ####
# Clean space
# rm(list=ls())

## Install and load packages
#install.packages("rstudioapi")
#install.packages("data.table")
library("rstudioapi")
library("data.table")
library("tidyr")

# Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "/01_scripts", replacement = "", x = current.path)
setwd(current.path)
rm(current.path)

## Info
# sessionInfo()

# Set variables
input.FN             <- "04_impute/all_inds_wgrs_and_panel_biallele_ai2.txt"
chr_indicator.string <- "NC_047"
output_folder        <- "04_impute/"

#### 01. Import data ####
input.dat <- fread(file = input.FN, sep = "\t", header = T)
input.dat[1:5,1:5]
dim(input.dat)

#### 02. Identify the unique chr present ####
## Identify the unique chr present based on chr_indicator string
chr.vec <- gsub(pattern = "\ .*", replacement = "", x = input.dat$mname)
chr.vec <- unique(chr.vec[grep(pattern = chr_indicator.string, x = chr.vec)])
chr.vec

#### 03. Subset the matrix into individual chr ####
# Subset the matrix into these chr
temp.mat <- NULL; matrix.list <- list(); output.FN <- NULL
for(i in 1:length(chr.vec)){
  
  # Reporting
  #print(paste0("Working on chr ", chr.vec[i]))
  
  # Obtain the section of the matrix with the target chr
  
  temp.mat <- input.dat[grep(pattern = chr.vec[i], x = input.dat$mname), ]
  
  # Reporting
  print(paste0("Number records for ", chr.vec[i], ": ", nrow(temp.mat)))
  
  temp.mat <- as.data.frame(temp.mat)
  dim(temp.mat)
  temp.mat[1:5,1:5]
  
  temp.mat <- t(temp.mat)
  temp.mat <- as.data.frame(temp.mat)
  dim(temp.mat)
  temp.mat[1:5,1:5]
  
  # Save to list
  matrix.list[[chr.vec[i]]] <- temp.mat
  
  # Drop mname row
  temp.mat <- temp.mat[grep(pattern = "mname", x = rownames(temp.mat), invert = T), ]
  temp.mat[1:5,1:5]
  
  # Write out
  output.FN <- paste0(output_folder, "ai2_input_", chr.vec[i], ".txt")
  #print(paste0("Writing out as ", output.FN))
  fwrite(x = temp.mat, file = output.FN, sep = " "
         , quote = F, col.names = F, row.names = T
         )
  
}

# Next: run ai2 on each of the output txt files
