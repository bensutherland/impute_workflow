# Plot GEMMA output
#  initialized 2024-07-26
#  Ben J. G. Sutherland (VIU), incl. code dev by Konstantin Divilov

#### 00. Front Matter ####
# Clear space
# rm(list=ls())

## Install packages
#devtools::install_github('kaustubhad/fastman',build_vignettes = TRUE)
#install.packages("missMethods")
#install.packages("tidyr")
#install.packages("ggplot2")
#install.packages("data.table)

## Load libraries
library(fastman)
library(missMethods)
library(tidyr)
library(ggplot2)
library(data.table)

## Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path) # take main directory
setwd(current.path)
rm(current.path)

# User set variables
gemma_output.FN       <- "07_GWAS/output/gwas_all_fam_covariate.assoc.txt"
phenotype_of_interest <- "survival_state" # DPE or survival_state
highlight_snps.FN <- "07_GWAS/denovo_snp_ids.txt"

# Read in data
gemma_gwas <- fread(file = gemma_output.FN, header = T)
gemma_gwas <- as.data.frame(gemma_gwas)
head(gemma_gwas)

highlight_snps.df <- read.delim(file = highlight_snps.FN, header = F, sep = "\t")
head(highlight_snps.df)

# Convert mname into chr and pos
gemma_gwas <- separate(data = gemma_gwas, col = "rs", into = c("chr", "pos"), sep = "__", remove = F)
head(gemma_gwas)
gemma_gwas$pos <- as.numeric(gemma_gwas$pos)

# convert to linkage group (LG), based on https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_902806645.1/
# LG to Chr info can be found in Sup File in https://doi.org/10.1093/gigascience/giab020
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

pdf(file = paste0("07_GWAS/Manhattan_plot_", phenotype_of_interest, ".pdf"), width = 11, height = 7)
par(mfrow = c(1,1), mar = c(5,4,4,2) +0.1, mgp = c(3,1,0))
fastman(m = gemma_gwas, chr = "chr", bp = "pos", p = "p_wald"
        , genomewideline = -log10(0.05/nrow(gemma_gwas))
        #, suggestiveline = -log10(0.05)
        , cex = 0.7, cex.lab = 1, cex.axis = 1
        , ylim = c(0,10)
        )
dev.off()



head(gemma_gwas)
# colnames(gemma_gwas)[which(colnames(gemma_gwas)=="rs")] <- "SNP"
head(gemma_gwas)

#highlight_vec <- gemma_gwas[which(gemma_gwas$p_wald==min(gemma_gwas$p_wald)), "SNP"]
highlight.vec <- highlight_snps.df$V1
length(highlight.vec)
highlight.vec <- highlight.vec[highlight.vec %in% gemma_gwas$rs]
length(highlight.vec)

# gemma_gwas$cex <- gemma_gwas$SNP %in% highlight.vec
# gemma_gwas$cex[which(gemma_gwas$cex==TRUE)] <- 1.5
# gemma_gwas$cex[which(gemma_gwas$cex==0)] <- 0.7
# table(gemma_gwas$cex)

pdf(file = paste0("07_GWAS/Manhattan_plot_", phenotype_of_interest, "_annotated.pdf"), width = 10.8, height = 5.2)
par(mfrow = c(1,1), mar = c(5,4,4,2) +0.1, mgp = c(3,1,0))
fastman(m = gemma_gwas, chr = "chr", bp = "pos", p = "p_wald", snp = "rs"
        , genomewideline = -log10(0.05/nrow(gemma_gwas))
        , suggestiveline = FALSE
        #, cex = 1.5
        #, cex = gemma_gwas$cex
        #, cex = ifelse(gemma_gwas$SNP %in% highlight_vec, 5, 1)
        , cex = 0.7
        , cex.lab = 1
        , cex.axis = 1
        #, ylim = c(0,10)
        , col = "Set2"
        , highlight = highlight.vec
        , annotateHighlight = T
        #, annotateTop=FALSE
        , annotationCol = "black"
        , annotationAngle= 55
)
dev.off()

# Write out table with whether SNP is from panel
gemma_gwas$panel_snp <- gemma_gwas$rs %in% highlight.vec
write.table(x = gemma_gwas, file = gsub(pattern = "\\.txt", replacement = "_w_annot.txt", x =  gemma_output.FN)
            , sep = "\t", quote = F, row.names = F, col.names = T
            )

# Write out as R object so can be easily loaded back
save.image(file = gsub(pattern = "\\.txt", replacement = ".RData", x = gemma_output.FN))
