# Read output of bcftools stats to summarize
# B. Sutherland (2024-09-04)
# Set the input folder containing the output of BCF comparison from bcftools stats (01_scripts/run_bcftools_stats.sh)
# Note: requires that "GCTs.txt", *_GCsS.txt and *_psd.txt is present in set 'input_folder'
# Note: requires that steps have been taken in the README before this script

### Front Matter ####
# Clean space
# rm(list=ls())

## Install and load packages
#install.packages("rstudioapi")
library(rstudioapi)
library(tidyr)
library(ggpubr)

# Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "/01_scripts", replacement = "", x = current.path)
setwd(current.path)
rm(current.path)

# Set user variables
#input_folder <- "05_compare/panel_vs_wgrs/" 
#input_folder <- "05_compare/ai2_vs_empirical/"
input_folder <- "05_compare/fi3_vs_empirical/"

options(scipen=999999999)


# Set include string to identify the individuals that should be in the summary 
remove_inds    <- TRUE
include_string <- "ASY2"
# note: if working with imputation and individuals were used with HD data they should be excluded here

# Plotting options
#plot_type <- "include_PSD" # Was per-site discordance calculated? Opts: "include_PSD" or "no_PSD"
# deprecated

# Build full filenames
input_GCTs.FN <- paste0(input_folder, "GCTs.txt")
input_GCsS.FN <- paste0(input_folder, "GCsS.txt")
input_PSD.FN  <- paste0(input_folder, "PSD.txt")


#### 01. R-squared value per sample (GCsS) ####
# Load data
concord_by_sample.df  <- read.delim2(file = input_GCsS.FN, header = T, sep = "\t")
concord_by_sample.df[1:5,1:5]

# Clean up df colnames
colnames(concord_by_sample.df) <- gsub(pattern = "^X\\.[0-9]\\.|^X\\.[0-9][0-9]\\.|^X\\.\\.", replacement = "", x = colnames(concord_by_sample.df))
concord_by_sample.df[1:5,1:5]

# Ensure dosage.r.squared is numeric
concord_by_sample.df$dosage.r.squared <- as.numeric(x = concord_by_sample.df$dosage.r.squared)

# Calculate per-sample r-value from r-squared value
concord_by_sample.df$dosage.r <- sqrt(concord_by_sample.df$dosage.r.squared)

# Remove cases with all missing data
empty_ind <- concord_by_sample.df[rowSums(concord_by_sample.df[, grep(pattern = "matches", x = colnames(concord_by_sample.df))])==0, "sample"]
print(paste0("Removing ", empty_ind, " due to no records present"))
nrow(concord_by_sample.df) # before any removal
concord_by_sample.df <- concord_by_sample.df[!(concord_by_sample.df$sample %in% empty_ind), ]
nrow(concord_by_sample.df) # after removal

# Remove inds if needed
if(isTRUE(remove_inds)){
  
  # Retain excluded samples for downstream inspection
  concord_by_sample_excluded.df <-  concord_by_sample.df[grep(pattern = include_string, x = concord_by_sample.df$sample, invert = T), ]
  print(paste0("Number individuals removed from the concordance table: ", nrow(concord_by_sample_excluded.df)))
  
  # Only keep those individuals to be included
  concord_by_sample.df <- concord_by_sample.df[grep(pattern = include_string, x = concord_by_sample.df$sample), ]
  paste0("Number of individuals being assessed after excluding drop inds: ", nrow(concord_by_sample.df))
  
  
}else{
  
  print("Not excluding any individuals")
  
}

# Plot histogram (r-squared)
pdf(file = paste0(input_folder, "hist_concord_by_sample_rsquared.pdf"), width = 7.6, height = 4)
hist(x = concord_by_sample.df$dosage.r.squared, las = 1, main = ""
     , xlab = expression(Concordance ~ per ~ sample ~ (allelic ~ dosage ~ R^2) )
     , breaks = 20
     , xlim = c(0,1)
     )
dev.off()

# Summary
summary(concord_by_sample.df$dosage.r.squared)
sd(concord_by_sample.df$dosage.r.squared)

# Also provide summary for excluded samples
if(isTRUE(remove_inds)){
  
  # Summary for excluded
  print("These are the individuals have been removed")
  print(summary(concord_by_sample_excluded.df$dosage.r.squared))
  print(sd(concord_by_sample_excluded.df$dosage.r.squared))
  
}

# # Optional: identify samples with r2 < 0.5
# concord_by_sample.df[concord_by_sample.df$dosage.r.squared < 0.5, "sample"]

# Plot histogram (r-value)
pdf(file = paste0(input_folder, "hist_concord_by_sample_r.pdf"), width = 7.6, height = 4)
hist(x = concord_by_sample.df$dosage.r, las = 1, main = ""
     , xlab = expression(Concordance ~ per ~ sample ~ (allelic ~ dosage ~ R) )
     , breaks = 20
     , xlim = c(0,1)
)
dev.off()

# Plot histogram (r-value) with ggplot2
hist_plot_rval <- concord_by_sample.df %>%
                  ggplot( aes(x=dosage.r)) +
                  geom_histogram( bins=40, fill="darkgrey", color="#e9ecef") +
                  theme_bw()+
                  labs(x = "Per sample allelic dosage (r-value)") + 
                  xlim(0.4, 1)
  
hist_plot_rval

# Summary
summary(concord_by_sample.df$dosage.r)
sd(concord_by_sample.df$dosage.r)

# Also provide summary for excluded samples, if removing inds
if(isTRUE(remove_inds)){
  
  print("These are the individuals have been removed")
  print(summary(concord_by_sample_excluded.df$dosage.r))
  print(sd(concord_by_sample_excluded.df$dosage.r))

}


#### 02. Percent concordant per sample (GCTs) ####
# Load data
concord_table.df  <- read.delim2(file = input_GCTs.FN, header = T, sep = "\t")
concord_table.df[1:5,1:5]

# Clean up df colnames
colnames(concord_table.df) <- gsub(pattern = "^X\\.[0-9]\\.|^X\\.[0-9][0-9]\\.|^X\\.\\.", replacement = "", x = colnames(concord_table.df))
colnames(concord_table.df) <- gsub(pattern = "\\.\\.\\.\\.", replacement = "-", x = colnames(concord_table.df))
concord_table.df[1:5,1:10]

# Remove the indiv with all missing data, if there is one present, as identified above
print(paste0("Also removing any indiv with all missing data (if present), as above: ", empty_ind))
nrow(concord_table.df)
concord_table.df <- concord_table.df[!(concord_table.df$sample %in% empty_ind), ] 
nrow(concord_table.df)

# Detect concordant colnames
columns <- colnames(concord_table.df)
columns <- as.data.frame(columns)
columns <- separate(data = columns, col = "columns", into = c("vcf1", "vcf2"), sep = "-", remove = F)
columns <- columns$columns[columns$vcf1==columns$vcf2]
concord_cols <- columns[!is.na(columns)] # get rid of NAs
concord_cols <- concord_cols[grep(pattern = "missing", x = concord_cols, invert = T)] # get rid of missing
concord_cols

# Detect missing-value colnames
missing_cols <- colnames(concord_table.df)[grep(pattern = "missing", x = colnames(concord_table.df))]
missing_cols

# Detect discordant colnames
all_but_discord_cols <- c("GCTs", "sample", concord_cols, missing_cols)
discord_cols <- colnames(concord_table.df)[!(colnames(concord_table.df) %in% all_but_discord_cols)]
discord_cols

## Summarize concordance per ind
# Prepare df
concord_table_summary.df <- as.data.frame(concord_table.df$sample)
colnames(concord_table_summary.df) <- c("indiv")
head(concord_table_summary.df)

# Count concord, discord, and missing
concord_table_summary.df$concord_count <- rowSums(concord_table.df[, concord_cols])
concord_table_summary.df$discord_count <- rowSums(concord_table.df[, discord_cols])
concord_table_summary.df$missing_count <- rowSums(concord_table.df[, missing_cols])

# Reporting
head(concord_table_summary.df)

# Check (OK)
table(rowSums(concord_table_summary.df[,c("concord_count", "missing_count", "discord_count")]))
# note: should only be one value, all same, and the number of SNPs expected

# Calculate percentage of complete records that are concordant (exclude missing)
concord_table_summary.df$prop_corr <- concord_table_summary.df$concord_count / (concord_table_summary.df$concord_count + concord_table_summary.df$discord_count)

# Remove inds if needed
if(isTRUE(remove_inds)){
  
  # Create separate df for the keep inds, or exclude inds
  concord_table_summary_excluded.df <- concord_table_summary.df[grep(pattern = include_string, x = concord_table_summary.df$indiv, invert = T), ]
  
  concord_table_summary.df <- concord_table_summary.df[grep(pattern = include_string, x = concord_table_summary.df$indiv), ]
  paste0("Number of individuals being assessed after excluding drop inds: ", nrow(concord_table_summary.df))
  
  
}else {
  
  print("Not excluding any inds")
  
}


# Plot histogram of proportion correct (excluding missing)
pdf(file = paste0(input_folder, "hist_concord_by_sample_proportion_corr.pdf")
    , width = 7.6, height = 4)
hist(x = concord_table_summary.df$prop_corr, las = 1, main = ""
     , xlab = "Proportion concordant (%)"
     , breaks = 20
     , xlim = c(0,1)
)
dev.off()

# Summary for included
summary(concord_table_summary.df$prop_corr)
sd(concord_table_summary.df$prop_corr, na.rm = T)

# Plot histogram (concordance) with ggplot2
hist_plot_prop_concord <- concord_table_summary.df %>%
  ggplot( aes(x=prop_corr)) +
  geom_histogram( bins=40, fill="darkgrey", color="#e9ecef") +
  theme_bw()+
  labs(x = "Per sample proportion concordant") + 
  xlim(0.5, 1)

hist_plot_prop_concord

# Also summarize excluded, if excluding
if(isTRUE(remove_inds)){
  
    # Summarize excluded for completeness
    print("These are the values for the excluded individuals")
    print(summary(concord_table_summary_excluded.df$prop_corr))
    print(sd(concord_table_summary_excluded.df$prop_corr, na.rm = T))

}

# Plot proportion correct by missing
pdf(file = paste0(input_folder, "ppn_concord_by_missing_empirical_data.pdf")
    , width = 7.6, height = 4)
plot(x = concord_table_summary.df$missing_count, y = concord_table_summary.df$prop_corr
     , xlab = "Per sample missing data", ylab = "Proportion concordant (%)"
     , las = 1
     )
dev.off()

# ggplot option
scatter_missing_by_prop_concord <- ggplot(concord_table_summary.df, aes(prop_corr, missing_count)) + 
                                        geom_point() + 
                                        theme_bw() + 
                                        labs(x = "Per sample proportion concordant", y = "Per sample missing") + 
                                        xlim(0.5, 1)
  
scatter_missing_by_prop_concord


#### 03. Per-site discordance (PSD) ####
if(file.exists(x = input_PSD.FN)){
  
  # Load data
  psd_table.df  <- read.delim2(file = input_PSD.FN, header = T, sep = "\t")
  psd_table.df[1:5,]
  
  # Clean up column names
  colnames(psd_table.df) <- gsub(pattern = "^X\\.[0-9]\\.|^X\\.[0-9][0-9]\\.|^X\\.\\.", replacement = "", x = colnames(psd_table.df))
  colnames(psd_table.df) <- gsub(pattern = "\\.\\.\\.\\.", replacement = "-", x = colnames(psd_table.df))
  psd_table.df[1:5,]
  
  # Create unique locus identifier
  psd_table.df$locus_id <- paste0(psd_table.df$CHROM, "_", psd_table.df$POS)
  psd_table.df$percent.corr <- psd_table.df$Number.of.matches / (psd_table.df$Number.of.matches + psd_table.df$Number.of.mismatches)
  psd_table.df$num_typed <- psd_table.df$Number.of.matches + psd_table.df$Number.of.mismatches
  
  pdf(file = paste0(input_folder, "per_locus_percent_corr.pdf"), width = 7.6, height = 4)
  hist(psd_table.df$percent.corr, main = "", xlab = "Percent correct (%)", las = 1)
  dev.off()
  
  if(nrow(psd_table.df) < 10000){
    
    pdf(file = paste0(input_folder, "per_locus_percent_corr_by_num_typed.pdf"), width = 7.6, height = 4)
    plot(y = psd_table.df$percent.corr, x = psd_table.df$num_typed, las = 1
         , ylab = "Percent correct (%)", xlab = "Number loci typed"
    )
    dev.off()
    
    
  }else{
    
    print("Too many loci to plot, not plotting percent corr by number typed to not crash the program.")
    
  }
  
  # How many loci have lower than cutoff percent correct? 
  bad_loci <- psd_table.df[psd_table.df$percent.corr < 0.5, "locus_id"]
  length(bad_loci)   # 'bad' loci
  nrow(psd_table.df) # All loci
  write.table(x = bad_loci, file = paste0(input_folder, "bad_loci.txt"), col.names = F, row.names = F)
  
  # Summary
  summary(psd_table.df$percent.corr)
  sd(psd_table.df$percent.corr, na.rm = T)
  
  psd_table_clean.df <- psd_table.df[psd_table.df$percent.corr >= 0.5, ]
  nrow(psd_table_clean.df)
  summary(psd_table_clean.df$percent.corr)
  sd(psd_table_clean.df$percent.corr, na.rm = T)
  
  # ggplot option
  psd_plot <- psd_table.df %>% 
                ggplot( aes(x = percent.corr)) + 
                geom_histogram( bins = 20, fill = "darkgrey", color = "#e9ecef") + 
                theme_bw() + 
                labs(x = "Per locus proportion concordant")
                #xlim(0, 1.2) # Does not render well, loses last bar
  psd_plot
  
  pdf(file = paste0(input_folder, "per_locus_percent_corr_ggplot_hist.pdf"), width = 7.6, height = 4)
  print(psd_plot)
  dev.off()
  
  
  
}else{
  
  print(paste0("No PSD file is present in ", input_folder))
  
}


##### End matter ####
save.image(file = paste0(input_folder, "assess_bcftools_stats_output_for_plots.RData"))

# Next, can plot using script '01_scripts/plot_impute_vs_empirical_multipanel.R'

# End of assessment of bcftools stats output
