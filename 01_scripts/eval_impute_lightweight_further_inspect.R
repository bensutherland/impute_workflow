# Further inspect the outcome of the concordance comparison
# B. Sutherland (2024-08-07)

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

## Info
# sessionInfo()

# Set variables
#per_chr.FN  <- "13_impute_compare_no_novel_no_MERR/concord_eval_NC_047559.1_comparison.txt"
all_data.FN <- "13_impute_compare/all_loci_data_comparison.txt" 

#### 01. Load data ####
# Read in per chr data
per_chr.df <- fread(file = per_chr.FN, sep = "\t")
dim(per_chr.df)
per_chr.df <- as.data.frame(per_chr.df) # convert to df
per_chr.df[1:5,1:5]

hist(per_chr.df$prop.match, breaks = 10, xlab = "Proportion match (%)", las = 1)

# which are the top samples for matching, and which are the bottom? 
per_chr.df <- per_chr.df[with(per_chr.df, order(per_chr.df$prop.match, decreasing = T)), ]
head(per_chr.df, n = 10)
tail(per_chr.df, n = 10)

# Read in all data
all_data.df <- fread(file = all_data.FN, sep = "\t")
dim(all_data.df)
all_data.df <- as.data.frame(all_data.df) # convert to df
all_data.df[1:5,1:5]

# Convert back to chr and marker
all_data.df <- separate(data = all_data.df, col = "mname", into = c("chr", "pos"), sep = "__", remove = F)
all_data.df$pos <- as.numeric(all_data.df$pos)

all_data_subset.df <- all_data.df[all_data.df$chr=="NC_047559.1", ]
dim(all_data_subset.df)

# They should arrive matched
colnames(all_data_subset.df)[4:21]
tail(colnames(all_data_subset.df))

# Order by position in chr
all_data_subset.df <- all_data_subset.df[with(all_data_subset.df, order(all_data_subset.df$pos, decreasing = F)), ]
all_data_subset.df[1:5,1:5]
dim(all_data_subset.df)


### ASIDE, check correlation to be sure ####
empirical.vector <- all_data_subset.df$ASY2_114_R1_1_empirical
imputed.vector   <- all_data_subset.df$ASY2_114_R1_1_imputed

empirical.vector <- gsub(pattern = "9", replacement = NA, x = empirical.vector)
imputed.vector   <- gsub(pattern = "9", replacement = NA, x = imputed.vector)

empirical.vector <- as.numeric(empirical.vector)
imputed.vector   <- as.numeric(imputed.vector)

length(empirical.vector)
length(imputed.vector)

cor(x = empirical.vector, y = imputed.vector, method = "pearson", use = "pairwise.complete.obs")

# Remove missing data 
all_data_subset_plot.df <- all_data_subset.df[all_data_subset.df$ASY2_114_R1_1_empirical!=9, c("ASY2_114_R1_1_empirical", "ASY2_114_R1_1_imputed")]
dim(all_data_subset_plot.df)

pdf(file = "13_impute_compare_no_novel_no_MERR/plot_err_emp-imp_NC_047559.1_114_R1_1.pdf", width = 20, height = 7)
plot(all_data_subset_plot.df$ASY2_114_R1_1_empirical - all_data_subset_plot.df$ASY2_114_R1_1_imputed
     )
dev.off()
