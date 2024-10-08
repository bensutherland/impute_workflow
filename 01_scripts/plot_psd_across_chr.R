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
per_locus_af.FN <- "05_compare_all_loci/fi3_vs_empirical/all_inds_empirical_shared_AF.txt" 
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


## Read in per-locus (empirical) AF
per_locus_af.df <- fread(file = per_locus_af.FN, sep = " ", quote = F, header = F)
per_locus_af.df <- as.data.frame(per_locus_af.df)
colnames(per_locus_af.df) <- c("chr", "pos", "AF")
per_locus_af.df$pos <- as.numeric(per_locus_af.df$pos)

# make ID
per_locus_af.df$locus.id <- paste0(per_locus_af.df$chr, "__", per_locus_af.df$pos)
head(per_locus_af.df)


#### MAF work ####
# MAF inspect
psd_and_af.df <- merge(x = psd.df, y = per_locus_af.df, by = "locus.id")
nrow(psd.df)
nrow(psd_and_af.df)
head(psd_and_af.df)

# Bin by AF
# Convert to MAF
# for(i in 1:nrow(psd_and_af.df)){
#   
#   if(psd_and_af.df$AF[i] > 0.5){
#     
#     psd_and_af.df$AF[i] <- 1 - psd_and_af.df$AF[i]
#     
#   }
#   
# }

MAF_bins <- seq(from = 0, to = 1.05, by = 0.05)

psd_and_af.df$bin <- NA
for(i in 2:length(MAF_bins)){
  
  psd_and_af.df[psd_and_af.df$AF < MAF_bins[i] & psd_and_af.df$AF > MAF_bins[i-1], "bin"] <- MAF_bins[i]
  
  psd_and_af.df[psd_and_af.df$AF == 1.0, "bin"] <- 1.0
  
}

table(psd_and_af.df$bin)
head(psd_and_af.df)

# Drop AF = 1.0
psd_and_af_filt.df <- psd_and_af.df[psd_and_af.df$AF!=1, ] 
dim(psd_and_af.df)
dim(psd_and_af_filt.df)

# With transparency (right)
p2 <- ggplot(data=psd_and_af_filt.df, aes(x=prop.match, group=bin, fill=bin)) +
  geom_density(adjust=1.5, alpha=.4) 
p2

pdf(file = "05_compare_all_loci/fi3_vs_empirical/boxplot_perloc_prop_match_by_AF.pdf", width = 12, height = 5)
boxplot(psd_and_af_filt.df$prop.match ~ psd_and_af_filt.df$bin, las = 1
        , xlab = "Allele frequency bin", ylab = "Per-locus proportion match"
        )
dev.off()


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
