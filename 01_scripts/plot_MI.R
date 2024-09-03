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
datatype <- "HD"
#datatype <- "LD"

input_MI_freq.FN <- paste0("06_screen_loci/", datatype, "_MI_freq.txt")

# Set the number of trios in total present
## TODO: naming variable
#MI_trios.FN <- "06_screen_loci/trios_limited.txt"
MI_trios.FN <- "06_screen_loci/trios.txt"

# Read in data
MI_freq.df <- read.table(file = input_MI_freq.FN, header = F, sep = " ")
colnames(MI_freq.df) <- c("chr", "pos", "num_trios_MI")
head(MI_freq.df)
num_loci <- nrow(MI_freq.df)

MI_trios.df <- read.table(file = MI_trios.FN, header = F, sep = ",")
num_trios <- nrow(MI_trios.df)

# Number of loci that have MI in at least half the trios
table(MI_freq.df$num_trios_MI >= (num_trios/2)) #31 of 5800937
table(MI_freq.df$num_trios_MI >= (num_trios/10)) #348329 of 5800937


output_histogram.FN <- paste0("06_screen_loci/histogram_of_MI_", datatype, ".pdf")

pdf(file = output_histogram.FN, width = 5, height = 4)
hist(MI_freq.df$num_trios_MI, las = 1, main = "", xlab = "Number trios displaying Mendelian incompatibility"
     , breaks = 15
     )

text(x = 110, y = 1200, labels = paste0("# loci > 100 MI: ", nrow(MI_freq.df[MI_freq.df$num_trios_MI > 100, ])))
dev.off()
     

# End
