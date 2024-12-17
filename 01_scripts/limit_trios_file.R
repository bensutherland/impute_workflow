# Limit trios file
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
trios.FN     <- "06_screen_loci/trios.txt"
present_samples.FN <- "06_screen_loci/wgrs_present_offspring.txt"


trios.df <- read.table(file = trios.FN, sep = ",")
head(trios.df)

present_samples <- read.table(file = present_samples.FN, header = F)
head(present_samples, n = 10)

trios.df <- trios.df[trios.df$V3 %in% present_samples$V1, ]
nrow(trios.df)
write.table(x = trios.df, file = "06_screen_loci/trios_limited.txt", col.names = F, row.names = F, sep = ","
            , quote = F)
