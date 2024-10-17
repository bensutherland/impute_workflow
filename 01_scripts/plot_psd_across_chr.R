# Plot per-site discordance (PSD) across chromosomes based on output of bcftools stats
#  initialized 2024-10-08
#  Ben J. G. Sutherland (VIU)

#### 00. Front Matter ####
# Clear space
# rm(list=ls())

## Install packages
#install.packages("tidyr")
#install.packages("ggplot2")
#install.packages("data.table)

## Load libraries
library(tidyr)
library(ggplot2)
library(data.table)

## Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path) # take main directory
setwd(current.path)
rm(current.path)

# User set variables
psd.FN            <- "05_compare_all_loci/fi3_vs_empirical/PSD.txt"
highlight_snps.FN <- "07_GWAS/denovo_snp_ids.txt"
# per_locus_af.FN <- "05_compare_all_loci/fi3_vs_empirical/all_inds_empirical_shared_AF.txt"
per_locus_maf.FN <- "04_impute_all_loci/all_inds_wgrs_and_panel_biallele_only_fi3_imputed_with_tags_perloc_MAF.txt" # imputed
#per_locus_maf.FN <- "05_compare_all_loci/panel_vs_wgrs/parents_and_offspring_wgrs_renamed_with_tags_perloc_MAF.txt"   # empirical
bin_size <- 100000
options(scipen=999999)


#### Read in data ####
## Read in per-site discordance data
psd.df <- fread(file = psd.FN)
psd.df <- as.data.frame(psd.df)
head(psd.df)
colnames(psd.df) <- c("psd", "chr", "pos", "num.match", "num.mismatch", "nrd")
head(psd.df)
str(psd.df)

# Calc prop match
psd.df$prop.match <- psd.df$num.match / (psd.df$num.mismatch + psd.df$num.match)

# Create ID
psd.df$locus.id <- paste0(psd.df$chr, "__", psd.df$pos)
head(psd.df)

## Read in chr and pos of de novo SNPs used (panel loci)
denovo_snps.df <- read.delim(file = highlight_snps.FN, header = F)
denovo_snps.df <- as.data.frame(denovo_snps.df)
denovo_snps.df <- separate(data = denovo_snps.df, col = "V1", into = c("chr", "pos"), sep = "__", remove = F)
denovo_snps.df$pos <- as.numeric(denovo_snps.df$pos)
head(denovo_snps.df)

# ## Read in per-locus (empirical) AF
# per_locus_af.df <- fread(file = per_locus_af.FN, sep = " ", quote = F, header = F)
# per_locus_af.df <- as.data.frame(per_locus_af.df)
# colnames(per_locus_af.df) <- c("chr", "pos", "AF")
# per_locus_af.df$pos <- as.numeric(per_locus_af.df$pos)
# 
# # make ID
# per_locus_af.df$locus.id <- paste0(per_locus_af.df$chr, "__", per_locus_af.df$pos)
# head(per_locus_af.df)

## Read in per-locus (imputed) MAF
per_locus_maf.df <- fread(file = per_locus_maf.FN, sep = " ", quote = F, header = F)
per_locus_maf.df <- as.data.frame(per_locus_maf.df)
colnames(per_locus_maf.df) <- c("chr", "pos", "MAF")
per_locus_maf.df$pos <- as.numeric(per_locus_maf.df$pos) # make numeric
per_locus_maf.df$locus.id <- paste0(per_locus_maf.df$chr, "__", per_locus_maf.df$pos) # make locus identifier
head(per_locus_maf.df)


#### MAF work ####
# MAF inspect
psd_and_maf.df <- merge(x = psd.df, y = per_locus_maf.df, by = "locus.id")
nrow(psd.df)
nrow(psd_and_maf.df)
head(psd_and_maf.df)


##### Convert AF to MAF #####
# Note: no longer needed, as generate MAF with bcftools
# # Make new column for MAF
# psd_and_af.df$MAF <- NA
# 
# # For each locus, check if greater than 0.5, and if so, subtract from 1 to make MAF
# for(i in 1:nrow(psd_and_af.df)){
# 
#   if(psd_and_af.df$AF[i] > 0.5){
# 
#     psd_and_af.df$MAF[i] <- 1 - psd_and_af.df$AF[i]
# 
# # Or else just keep the AF as MAF if less than or equal to 0.5
#   }else if(psd_and_af.df$AF[i] <= 0.5){
#     
#     psd_and_af.df$MAF[i] <- psd_and_af.df$AF[i]
#     
#   }
# 
# }


# Create bins for AF 
MAF_bins <- seq(from = 0, to = 0.55, by = 0.05)

# Designate the bin per locus
psd_and_maf.df$bin <- NA

# If MAF == 0, it is monomorphic, bin it in the 0 bin
psd_and_maf.df[psd_and_maf.df$MAF == 0, "bin"] <- 0

for(i in 2:length(MAF_bins)){
  
  print(paste0("Adding bin attribute for MAF bin ", MAF_bins[i]))
  
  psd_and_maf.df[psd_and_maf.df$MAF > MAF_bins[i-1] & psd_and_maf.df$MAF <= MAF_bins[i], "bin"] <- MAF_bins[i]
  
}

# Summarize
table(psd_and_maf.df$bin, useNA = "ifany")
head(psd_and_maf.df)
#head(psd_and_maf.df[is.na(psd_and_maf.df$bin), ] )

# Drop monomorphs
#psd_and_maf_filt.df <- psd_and_maf.df[psd_and_maf.df$MAF > 0, ] 
psd_and_maf_filt.df <- psd_and_maf.df
dim(psd_and_maf_filt.df)

# With transparency (right)
p2 <- ggplot(data=psd_and_maf_filt.df, aes(x=prop.match, group=bin, fill=bin)) +
  geom_density(adjust=1.5, alpha=.4) 
p2

p2 <- ggplot(data=psd_and_maf_filt.df, aes(x=prop.match, group=bin, fill=bin)) +
  geom_density(adjust=1.5, alpha=.4) +
  facet_wrap(~bin)
p2

aggregate(prop.match ~ bin , psd_and_maf_filt.df, mean)


pdf(file = "05_compare_all_loci/fi3_vs_empirical/boxplot_perloc_prop_match_by_AF.pdf", width = 12, height = 5)
boxplot(psd_and_af_filt.df$prop.match ~ psd_and_af_filt.df$bin, las = 1
        , xlab = "Allele frequency bin", ylab = "Per-locus proportion match"
        )
dev.off()


# What would happen if you filtered under a certain MAF?
# currently use < 0.3 for imputed, <= 0.15 for empirical
test.df <- psd_and_maf.df[psd_and_maf.df$MAF > 0 & psd_and_maf.df$MAF < 0.3 , ] 
dim(test.df)
summary(test.df$prop.match)
sd(test.df$prop.match)

# Compare to all data, no monomorphs
summary(psd_and_maf.df[psd_and_maf.df$MAF > 0 , "prop.match"]) 
sd(psd_and_maf.df[psd_and_maf.df$MAF > 0 , "prop.match"]) 
length(psd_and_maf.df[psd_and_maf.df$MAF > 0 , "prop.match"]) # number records


#### Per chromosome plotting ####
present_chr <- unique(psd.df$chr)
present_chr

chr.oi <- NULL; psd_subset.df <- NULL; panel_subset.df <- NULL
max_pos <- NULL; sequence_of_bins <- NULL; psd_summary.df <- NULL
for(c in 1:length(present_chr)){
  
  # Select the chr
  chr.oi <- present_chr[c]
  print(paste0("Working on chr: ", chr.oi))
  
  # Subset to the chr of interest
  psd_subset.df <- psd.df[psd.df$chr == chr.oi, ]
  
  # Reporting
  print(paste0("This chr has ", nrow(psd_subset.df), " loci"))
  
  # Keep the panel info for this chr
  panel_subset.df <- denovo_snps.df[denovo_snps.df$chr==chr.oi, ] 
  
  # Reporting
  print(paste0("This chr has ", nrow(panel_subset.df), " loci"))
  
  # Reporting
  print("Creating positional bins across chr")
  
  # Find the max position (plus buffer)
  max_pos <- max(psd.df$pos) + bin_size
  
  # Determine what bins will be included
  sequence_of_bins <- seq(from = 1, to = max_pos, by = bin_size)
  
  # Provide vector for bin
  psd_subset.df$bin <- NA
  
  # Assign each locus to appropriate bin
  print("Assigning loci to appropriate bins")
  
  for(i in 2:length(sequence_of_bins)){
    
    # For those loci within this bin size (less than i and greater than i-1), assign the bin name
    psd_subset.df[psd_subset.df$pos < sequence_of_bins[i] & psd_subset.df$pos > sequence_of_bins[i-1], "bin"] <- sequence_of_bins[i]
    
  }
  
  # Calculate the average of the prop match for each bin
  psd_summary.df <- aggregate(prop.match ~ bin, data = psd_subset.df, FUN = mean)
  
  # Plot output
  pdf(file = paste0("05_compare/avg_prop_corr_across_chr_", chr.oi, ".pdf"), width = 12, height = 5)
  plot(psd_summary.df$prop.match ~ psd_summary.df$bin, type = "l", lty  = 1, las = 1
       , ylab = "Avg. per locus prop. match"
       , xlab = "Position (bp)"
       , ylim = c(0.5,1)
       , main = paste0("chr: ", chr.oi, "; ", nrow(psd_subset.df), " loci; ", nrow(panel_subset.df), " panel loci")
  )
  
  
  # Add vertical separators
  abline(v = panel_subset.df$pos, lty = 3, col = "grey60")
  dev.off()
  
}




# End
