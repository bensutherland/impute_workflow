# Read in Fimpute3 file and convert to ai2 file
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
fi3_input.FN     <- "../FImpute3_test/output/genotypes_imp.txt"

fi3.df <- fread(file =  fi3_input.FN, sep = "\t")
fi3.df <- as.data.frame(fi3.df)
dim(fi3.df)
#str(fi3.df)

fi3_markers.FN     <- "../FImpute3_test/output/snp_info.txt" 

markers.df <- fread(file = fi3_markers.FN, sep = "\t")
markers.df <- as.data.frame(markers.df)
dim(markers.df)
head(markers.df)

genos_imputed.df <- fi3.df$ID
genos_imputed.df <- as.data.frame(genos_imputed.df)
colnames(genos_imputed.df) <- "ind"
head(genos_imputed.df)

output.mat <- matrix(data = NA, nrow = nrow(genos_imputed.df), ncol = nchar(fi3.df$Calls...[1]))

nucl.vec <- NULL; storage.list  <- list();
for(i in 1:nrow(fi3.df)){
  
   # Separate the genotypes into individual characters
   nucl.vec <- unlist(strsplit(x = fi3.df$Calls...[i], split = ""))
   
   # Add this into the collector matrix
   output.mat[i,] <- nucl.vec
  
}

output.df <- as.data.frame(output.mat)
dim(output.df)

colnames(x = output.df) <- markers.df$SNPID
output.df[1:4,1:5]

rownames(x = output.df) <- genos_imputed.df$ind
output.df[1:4,1:5]

fwrite(x = output.df, file = "05_compare/fi3_output_converted.txt", sep = "\t", row.names = T)



