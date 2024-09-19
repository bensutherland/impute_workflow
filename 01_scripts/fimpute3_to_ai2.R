# Read in chr-separated Fimpute3 output files and convert to ai2 format
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
fi3_folder  <- "04_impute/fimpute/" # the folder that holds the automated output directories for fimpute3
imputed.FN  <- "genotypes_imp.txt"  # automatic naming of the imputed genotypes by fi3
snp_info.FN <- "snp_info.txt"       # automatic naming of the SNP info (map) file

# What are the output folders? 
output.dirs <- list.dirs(path = fi3_folder, full.names = F)
output.dirs <- output.dirs[grep(pattern = "output", x = output.dirs)]

# output.dirs[1]
# list.files(path = paste0(fi3_folder, output.dirs[1]))

# Read in Fi3 imputed output genotypes and map file within an output folder
chr.num <- NULL; per_chr_data.list <- list(); imputed_genos.df <- NULL; markers.df <- NULL; all_results.df <- NULL
output.mat <- NULL; output.FN <- NULL
for(i in 1:length(output.dirs)){
  
  # Identify which numbered chr this round, based on the directory going in
  chr.num <- gsub(pattern = "output_folder_", replacement = "", x =  output.dirs[i])
  
  # Create output filename
  output.FN <- paste0(fi3_folder, "fi3_converted_chr_", chr.num, ".txt")
  
  # Read in imputed genos
  imputed_genos.df <- fread(file = paste0(fi3_folder, output.dirs[i], "/", imputed.FN), sep = "\t")
  imputed_genos.df <- as.data.frame(imputed_genos.df)
  dim(imputed_genos.df)
  str(imputed_genos.df)
  
  # Read in snp info
  markers.df <- fread(file = paste0(fi3_folder, output.dirs[i], "/", snp_info.FN), sep = "\t")
  markers.df <- as.data.frame(markers.df)
  dim(markers.df)
  head(markers.df)
  
  # Start an all results collector file
  all_results.df <- imputed_genos.df$ID
  all_results.df <- as.data.frame(all_results.df)
  colnames(all_results.df) <- "ind"
  
  # Create an empty matrix to be filled by genotypes
  output.mat <- matrix(data = NA, nrow = nrow(all_results.df), ncol = nchar(imputed_genos.df$Calls...[1])) # note: can be any of the rows' calls (hence the 1), they should all be the same
  
  # Loop over to pull out the genotypes from the single, no-space sep string
  nucl.vec <- NULL; storage.list  <- list()
  for(j in 1:nrow(imputed_genos.df)){
    
    # Separate the genotypes into individual characters
    nucl.vec <- unlist(strsplit(x = imputed_genos.df$Calls...[j], split = ""))
    
    # Add this into the collector matrix
    output.mat[j,] <- nucl.vec
    
  }
  
  # Make the collector matrix a df
  output.df <- as.data.frame(output.mat)
  dim(output.df)
  
  # And add back in the marker names
  colnames(x = output.df) <- markers.df$SNPID
  output.df[1:4,1:5]
  
  # Add back in the sample names
  rownames(x = output.df) <- all_results.df$ind
  output.df[1:4,1:5]
 
  
  
  output_t.df <- t(output.df)
  dim(output_t.df)
  output_t.df <- as.data.frame(output_t.df)
  output_t.df[1:5,1:5]
  
  fwrite(x = output_t.df, file = output.FN, sep = "\t", row.names = T)
  
  per_chr_data.list[[i]] <- output_t.df 
  
}


# Need to join all the chr back together
all_chr.df <- bind_rows(per_chr_data.list
                     #, .id = "column_label"
                     )
dim(all_chr.df)
all_chr.df[1:5,1:5]

# Bring rownames into mname
all_chr.df$mname <- rownames(all_chr.df)

# Put mname as first column, followed by everything else
all_chr.df <- all_chr.df %>% select(mname, everything())
all_chr.df[1:5,1:5]

# Remove double-underscore separator to make it more similar to ai2 output
all_chr.df$mname <- gsub(pattern = "__", replacement = " ", x = all_chr.df$mname)

# Set output filename
all_output.FN <- paste0(fi3_folder, "fi3_loci_by_inds_all_imputed_chr.txt")

# Write out full dataset, in ai2-format
fwrite(x = all_chr.df, file = all_output.FN, sep = "\t", row.names = F)
