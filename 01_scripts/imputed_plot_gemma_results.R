# Analyze gemma files
#  requires (...)
#  initialized 2024-07-26
#  Ben J. G. Sutherland (VIU), incl. code dev by Konstantin Divilov

#### 00. Front Matter ####
# Clear space
# rm(list=ls())

# Load libraries
#devtools::install_github('kaustubhad/fastman',build_vignettes = TRUE)
library(fastman)
#install.packages("missMethods")
library(missMethods)
library(tidyr)
library(ggplot2)

## Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path) # take main directory
setwd(current.path)
rm(current.path)

# User set variables
gemma_output.FN       <- "14_imputed_gwas/output/gwas_imputed_covar.assoc.txt"

gemma_gwas <- fread(file = gemma_output.FN, header = T)
head(gemma_gwas)

gemma_gwas <- as.data.frame(gemma_gwas)
gemma_gwas <- separate(data = gemma_gwas, col = "rs", into = c("chr.true", "pos.true"), sep = "__", remove = F)
head(gemma_gwas)
gemma_gwas$pos.true <- as.numeric(gemma_gwas$pos.true)

# convert to linkage group (LG), based on https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_902806645.1/
# LG to Chr info can be found in Sup File in https://doi.org/10.1093/gigascience/giab020
gemma_gwas$chr <- gemma_gwas$chr.true
gemma_gwas$chr <- gsub(pattern = "NC_047559.1", replacement = "Chr07", x = gemma_gwas$chr)
gemma_gwas$chr <- gsub(pattern = "NC_047560.1", replacement = "Chr01", x = gemma_gwas$chr)
gemma_gwas$chr <- gsub(pattern = "NC_047561.1", replacement = "Chr09", x = gemma_gwas$chr)
gemma_gwas$chr <- gsub(pattern = "NC_047562.1", replacement = "Chr06", x = gemma_gwas$chr)
gemma_gwas$chr <- gsub(pattern = "NC_047563.1", replacement = "Chr03", x = gemma_gwas$chr)
gemma_gwas$chr <- gsub(pattern = "NC_047564.1", replacement = "Chr02", x = gemma_gwas$chr)
gemma_gwas$chr <- gsub(pattern = "NC_047565.1", replacement = "Chr04", x = gemma_gwas$chr)
gemma_gwas$chr <- gsub(pattern = "NC_047566.1", replacement = "Chr05", x = gemma_gwas$chr)
gemma_gwas$chr <- gsub(pattern = "NC_047567.1", replacement = "Chr10", x = gemma_gwas$chr)
gemma_gwas$chr <- gsub(pattern = "NC_047568.1", replacement = "Chr08", x = gemma_gwas$chr)

gemma_gwas <- gemma_gwas[with(gemma_gwas, order(gemma_gwas$chr)), ]

hist(gemma_gwas$p_wald, breaks = 20)


par(mfrow = c(1,1), mar = c(5,4,4,2) +0.1, mgp = c(3,1,0))
fastman(m = gemma_gwas, chr = "chr", bp = "pos.true", p = "p_wald"
        , genomewideline = -log10(0.05/nrow(gemma_gwas))
        #, suggestiveline = -log10(0.05)
        , cex = 0.7, cex.lab = 1, cex.axis = 1
        , ylim = c(0,10)
        )




