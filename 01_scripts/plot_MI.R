# Plot distribution of MI
# B. Sutherland (2024-09-03)

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
MI_freq.FN     <- "06_screen_loci/MI_freq.txt" 

# Read in data
MI_freq.df <- read.table(file = MI_freq.FN, header = F, sep = " ")
colnames(MI_freq.df) <- c("chr", "pos", "num_trios_MI")
head(MI_freq.df)

pdf(file = "06_screen_loci/histogram_of_MI.pdf", width = 5, height = 4)
hist(MI_freq.df$num_trios_MI, las = 1, main = "", xlab = "Number trios displaying Mendelian incompatibility"
     , breaks = 15
     )

text(x = 110, y = 1200, labels = paste0("# loci > 100 MI: ", nrow(MI_freq.df[MI_freq.df$num_trios_MI > 100, ])))
dev.off()
     

# End
