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

# Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "/01_scripts", replacement = "", x = current.path)
setwd(current.path)
rm(current.path)

# Set user variables
#input_folder <- "05_compare/panel_vs_wgrs/" 
input_folder <- "05_compare/ai2_vs_empirical/"
#input_folder <- "05_compare/fi3_vs_empirical/"


# Set include string to exclude any individuals that should not be in the summary 
# ...for example, if individuals were used with HD data, and still are present here, exclude them
include_string <- "ASY2"

# TODO: optional - include list of samples to drop, if can't exclude using include_string


# Build full filenames
input_GCTs.FN <- paste0(input_folder, "GCTs.txt")
input_GCsS.FN <- paste0(input_folder, "GCsS.txt")
input_PSD.FN  <- paste0(input_folder, "PSD.txt")

# Not yet implemented
#r_type <- "not_squared"


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

# Only keep individuals with include_string, to exclude those that should not be in calculations
# ..but first, retain excluded samples for inspecting
concord_by_sample_excluded.df <-  concord_by_sample.df[grep(pattern = include_string, x = concord_by_sample.df$sample, invert = T), ]
print(paste0("Number individuals removed from the concordance table: ", nrow(concord_by_sample_excluded.df)))

# Only keep those individuals to be included
concord_by_sample.df <- concord_by_sample.df[grep(pattern = include_string, x = concord_by_sample.df$sample), ]
paste0("Number of individuals being assessed after excluding drop inds: ", nrow(concord_by_sample.df))

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

# Summary for excluded
print("These are the individuals have been removed")
summary(concord_by_sample_excluded.df$dosage.r.squared)
sd(concord_by_sample_excluded.df$dosage.r.squared)


# # Samples with r2 < 0.5
# concord_by_sample.df[concord_by_sample.df$dosage.r.squared < 0.5, "sample"]

# Plot histogram (r-value)
pdf(file = paste0(input_folder, "hist_concord_by_sample_r.pdf"), width = 7.6, height = 4)
hist(x = concord_by_sample.df$dosage.r, las = 1, main = ""
     , xlab = expression(Concordance ~ per ~ sample ~ (allelic ~ dosage ~ R) )
     , breaks = 20
     , xlim = c(0,1)
)
dev.off()

# Summary
summary(concord_by_sample.df$dosage.r)
sd(concord_by_sample.df$dosage.r)

# Summary for excluded
print("These are the individuals have been removed")
summary(concord_by_sample_excluded.df$dosage.r)
sd(concord_by_sample_excluded.df$dosage.r)


#### 02. Percent concordant per sample (GCTs) ####
# Load data
concord_table.df  <- read.delim2(file = input_GCTs.FN, header = T, sep = "\t")
concord_table.df[1:5,1:5]

# Clean up df colnames
colnames(concord_table.df) <- gsub(pattern = "^X\\.[0-9]\\.|^X\\.[0-9][0-9]\\.|^X\\.\\.", replacement = "", x = colnames(concord_table.df))
colnames(concord_table.df) <- gsub(pattern = "\\.\\.\\.\\.", replacement = "-", x = colnames(concord_table.df))
concord_table.df[1:5,1:10]

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

# Create separate df for the keep inds, or exclude inds
concord_table_summary_excluded.df <- concord_table_summary.df[grep(pattern = include_string, x = concord_table_summary.df$indiv, invert = T), ]

concord_table_summary.df <- concord_table_summary.df[grep(pattern = include_string, x = concord_table_summary.df$indiv), ]
paste0("Number of individuals being assessed after excluding drop inds: ", nrow(concord_table_summary.df))



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

# Summarize excluded for completeness
print("These are the values for the excluded individuals")
summary(concord_table_summary_excluded.df$prop_corr)
sd(concord_table_summary_excluded.df$prop_corr, na.rm = T)

# Plot proportion correct by missing
pdf(file = paste0(input_folder, "ppn_concord_by_missing_empirical_data.pdf")
    , width = 7.6, height = 4)
plot(x = concord_table_summary.df$missing_count, y = concord_table_summary.df$prop_corr
     , xlab = "Per sample missing data", ylab = "Proportion concordant (%)"
     , las = 1
     )
dev.off()


#### 03. Per-site discordance (PSD) ####
if(exists(x = input_PSD.FN)){
  
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
  
  pdf(file = "05_compare/panel_vs_wgrs/per_locus_percent_corr.pdf", width = 7.6, height = 4)
  hist(psd_table.df$percent.corr, main = "", xlab = "Percent correct (%)", las = 1)
  dev.off()
  
  pdf(file = "05_compare/panel_vs_wgrs/per_locus_percent_corr_by_num_typed.pdf", width = 7.6, height = 4)
  plot(y = psd_table.df$percent.corr, x = psd_table.df$num_typed, las = 1
       , ylab = "Percent correct (%)", xlab = "Number loci typed"
  )
  dev.off()
  
  # How many loci have lower than cutoff percent correct? 
  bad_loci <- psd_table.df[psd_table.df$percent.corr < 0.5, "locus_id"]
  length(bad_loci)
  write.table(x = bad_loci, file = "05_compare/panel_vs_wgrs/bad_loci.txt", col.names = F, row.names = F)
  nrow(psd_table.df)
  
  # Summary
  summary(psd_table.df$percent.corr)
  sd(psd_table.df$percent.corr, na.rm = T)
  
  psd_table_clean.df <- psd_table.df[psd_table.df$percent.corr >= 0.5, ]
  nrow(psd_table_clean.df)
  summary(psd_table_clean.df$percent.corr)
  sd(psd_table_clean.df$percent.corr, na.rm = T)
  
}else{
  
  print(paste0("No PSD file is present in ", input_folder))
  
}

# End of assessment of bcftools stats output
