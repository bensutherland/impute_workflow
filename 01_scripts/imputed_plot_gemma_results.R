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
#gemma_output.FN       <- "07_GWAS/ai2_imputed_survival_state_pheno/output/gwas_all_fam_covariate.assoc.txt"
gemma_output.FN       <- "07_GWAS/output/gwas_all_fam_covariate.assoc.txt"
phenotype_of_interest <- "survival_state" # DPE or survival_state
highlight_snps.FN <- "07_GWAS/denovo_snp_ids.txt"
caution_zones.FN <- "chr_locs_w_multimapper_loci.txt"  # This is if you are using multimappers to plot

plotting_multimappers <- FALSE # Set as TRUE if you want to plot multimappers and have the files available

plot_maxP <- 9999999 # set the maxP cutoff for the Fastman plot


#### 01. Read in data ####
# Read in GWAS output
gemma_gwas <- fread(file = gemma_output.FN, header = T)
gemma_gwas <- as.data.frame(gemma_gwas)
head(gemma_gwas)

# Read in the SNPs to highlight (single column with marker name, e.g., 'chr__position')
highlight_snps.df <- read.delim(file = highlight_snps.FN, header = F, sep = "\t")
head(highlight_snps.df)

# Read in multimapper caution zones if using
if(plotting_multimappers == TRUE){
  
  caution_zones.df <- read.delim(file = caution_zones.FN, header = T, sep = "\t")
  head(caution_zones.df)
  nrow(caution_zones.df)
  
}


#### 02. Format data as needed ####
# Convert gemma output file mname into chr and pos
gemma_gwas <- separate(data = gemma_gwas, col = "rs", into = c("chr", "pos"), sep = "__", remove = F)
head(gemma_gwas)
gemma_gwas$pos <- as.numeric(gemma_gwas$pos)

# Convert to linkage group (LG), based on https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_902806645.1/
# LG to chr info can be found in Sup File in https://doi.org/10.1093/gigascience/giab020
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


#### 03. Plotting ####
output.FN <- gsub(pattern = "\\/output.*", replacement = "", x = gemma_output.FN)

# Plot p-val histogram
pdf(file = paste0(output.FN, "/pval_hist_", phenotype_of_interest, ".pdf"), width = 7.5, height = 4.5)
hist(gemma_gwas$p_wald, breaks = 20)
dev.off()

# Plot Manhattan plot (default)
pdf(file = paste0(output.FN, "/Manhattan_plot_", phenotype_of_interest, ".pdf"), width = 9, height = 3.5)
par(mfrow = c(1,1), mar = c(5,4,4,2) +0.1, mgp = c(3,1,0))
fastman(m = gemma_gwas, chr = "chr", bp = "pos", p = "p_wald"
        , genomewideline = -log10(0.05/nrow(gemma_gwas))
        #, suggestiveline = -log10(0.05)
        , cex = 0.7, cex.lab = 1, cex.axis = 1
        #, ylim = c(0,10)
        , maxP = plot_maxP
        )
dev.off()


# Create highlight SNP vector, making sure that they are actually in the gemma output, and accounting
highlight.vec <- highlight_snps.df$V1
length(highlight.vec)
highlight.vec <- highlight.vec[highlight.vec %in% gemma_gwas$rs]
length(highlight.vec)


pdf(file = paste0(output.FN, "/Manhattan_plot_", phenotype_of_interest, "_annotated.pdf"), width = 10.8, height = 5.2)
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
        , maxP = plot_maxP
)
dev.off()



#### 04. Result retention ####
# Write out table with whether SNP is from panel
gemma_gwas$panel_snp <- gemma_gwas$rs %in% highlight.vec
write.table(x = gemma_gwas, file = gsub(pattern = "\\.txt", replacement = "_w_annot.txt", x =  gemma_output.FN)
            , sep = "\t", quote = F, row.names = F, col.names = T
            )


#### 05. Plotting multimapper caution zones ####
if(plotting_multimappers == TRUE){
  
  
        # Create identifier for matching
        caution_zones.vec <- paste0(caution_zones.df$chr, "__", caution_zones.df$pos)
        
        # How many directly are observed in the GEMMA output?
        caution_zones.vec[caution_zones.vec %in% gemma_gwas$rs]
        
        # Add some details to the caution zones df including up and downstream from the locus position
        head(caution_zones.df)
        str(caution_zones.df)
        caution_zones.df$upstream_loc <- caution_zones.df$pos - 100
        caution_zones.df$downstream_loc <- caution_zones.df$pos + 100
        head(caution_zones.df)
        
        # Separate the chr (scaffold) and position (scaff.pos) back in the GEMMA output
        gemma_gwas <- separate(data = gemma_gwas, col = "rs"
                               , into = c("scaffold", "scaff.pos"), sep = "__", remove = F
                               )
        gemma_gwas$scaff.pos <- as.numeric(gemma_gwas$scaff.pos) # ensure numeric
        
        # For each row of the caution zone df, identify if there are any GEMMA results in that section, and 
        #  if so, record their names
        head(gemma_gwas)
        gemma_gwas$caution_zone <- NA
        
        chr.oi <- NULL; slice <- NULL; problem_loci <- NULL; problem_vec <- NULL
        for(i in 1:nrow(caution_zones.df)){
          
          # Identify the chromosome of interest
          chr.oi <- caution_zones.df$chr[i]
          
          # Isolate the GEMMA results to only this chr
          slice <- gemma_gwas[gemma_gwas$scaffold==chr.oi, ]
          
          # Identify loci that are within the problem zones in this slice
          problem_loci <- slice[slice$scaff.pos >= caution_zones.df$upstream_loc[i] & slice$scaff.pos <= caution_zones.df$downstream_loc[i], "rs"] 
          
          # Retain all problem loci in a vector for highlighting below
          problem_vec <- c(problem_loci, problem_vec)
          
        }
        
        
        length(problem_vec)
        
        # Plot
        pdf(file = paste0("07_GWAS/Manhattan_plot_", phenotype_of_interest, "_annotated_with_multimapper_zones.pdf"), width = 10.8, height = 5.2)
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
                , highlight = problem_vec
                , annotateHighlight = T
                #, annotateTop=FALSE
                , annotationCol = "black"
                , annotationAngle= 55
                , maxP = plot_maxP
        )
        dev.off()
}



#### 06. Save out final results ####
# Write out as R object so can be easily loaded back
save.image(file = paste0(output.FN, "/post-plot.RData"))
# Downstream will use gemma_gwas object to multi-pannel plot GWAS output

