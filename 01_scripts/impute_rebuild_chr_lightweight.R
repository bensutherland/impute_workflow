# Rebuild the full datafile after ai2 that has been sep by chr
# B. Sutherland (2024-07-22)

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

# Set user variables
input_folder  <- "04_impute/"
output_folder <- "05_compare/" 

# Identify input files
input_files.vec <- list.files(path = input_folder, pattern = ".genotypes")

#### Read in imputed chr data ####
input.FN <- NULL; imputed_chr.df <- NULL
for(i in 1:length(input_files.vec)){
  
  # Set filename
  input.FN <- paste0(input_folder, input_files.vec[i])
  output.FN <- paste0(output_folder, basename(input.FN), "_transposed_to_combine.txt")
  
  # Reporting
  print(paste0("Reading in file ", input.FN))
  print(paste0("To write out as ", output.FN))
  
  # Read in chr data
  imputed_chr.df <- fread(file = input.FN)
  imputed_chr.df <- as.data.frame(x = imputed_chr.df)
  
  # Transpose to match ai2 input format
  imputed_chr.df <- t(imputed_chr.df)
  imputed_chr.df <- as.data.frame(x = imputed_chr.df)
  
  # If not the first chr, drop the ind row
  if(i==1){
    
    print("Retaining indnames for first chr.")
    
  }else if(i > 1){
    
    print("Dropping the ind name row")
    
    # Drop the ind row
    imputed_chr.df <- imputed_chr.df[-1,]
    
  }
  
  # Write out
  fwrite(x = imputed_chr.df, file = output.FN, sep = "\t", quote = F, col.names = F)
  
}

# Next, run "01_scripts/combine_transposed_ai2_output_and_mnames.sh"
