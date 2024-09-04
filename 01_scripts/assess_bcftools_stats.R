# Read output of bcftools stats to summarize
# B. Sutherland (2024-09-04)

### Front Matter ####
# Clean space
# rm(list=ls())

## Install and load packages
#install.packages("rstudioapi")
library(rstudioapi)

# Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "/01_scripts", replacement = "", x = current.path)
setwd(current.path)
rm(current.path)

# Set user variables
input_GCTs.FN <- "05_compare/panel_vs_wgrs/all_inds_wgrs_panel_comp_GCTs.txt"
input_GCsS.FN <- "05_compare/panel_vs_wgrs/all_inds_wgrs_panel_comp_GCsS.txt"

#### 01. R-squared value per sample ####
# Load data
concord_by_sample.df  <- read.delim2(file = input_GCsS.FN, header = T, sep = "\t")
concord_by_sample.df[1:5,1:5]

# Clean up column names
colnames(concord_by_sample.df) <- gsub(pattern = "^X\\.[0-9]\\.|^X\\.[0-9][0-9]\\.|^X\\.\\.", replacement = "", x = colnames(concord_by_sample.df))
concord_by_sample.df[1:5,1:5]

# Convert dosage rsquared to numeric
concord_by_sample.df$dosage.r.squared <- as.numeric(x = concord_by_sample.df$dosage.r.squared)

# Plot histogram
pdf(file = "05_compare/panel_vs_wgrs/hist_geno_concord_by_sample_r2.pdf", width = 7.6, height = 4)
hist(x = concord_by_sample.df$dosage.r.squared, las = 1, main = ""
     , xlab = "Genotype concordance by sample, dosage r2"
     , breaks = 20
     )
dev.off()


# Summary
length(concord_by_sample.df$sample)
summary(concord_by_sample.df$dosage.r.squared)
sd(concord_by_sample.df$dosage.r.squared)

# Samples with r2 < 0.5
concord_by_sample.df[concord_by_sample.df$dosage.r.squared < 0.5, "sample"]


#### 02. Percent concordant per sample ####
# Load data
concord_table.df  <- read.delim2(file = input_GCTs.FN, header = T, sep = "\t")
concord_table.df[1:5,1:5]

# Clean up column names
colnames(concord_table.df) <- gsub(pattern = "^X\\.[0-9]\\.|^X\\.[0-9][0-9]\\.|^X\\.\\.", replacement = "", x = colnames(concord_table.df))
colnames(concord_table.df) <- gsub(pattern = "\\.\\.\\.\\.", replacement = "-", x = colnames(concord_table.df))
concord_table.df[1:5,1:10]

# Detect concordant colnames
columns <- colnames(concord_table.df)
columns <- as.data.frame(columns)
columns <- separate(data = columns, col = "columns", into = c("vcf1", "vcf2"), sep = "-", remove = F)
columns <- columns$columns[columns$vcf1==columns$vcf2]
concord_cols <- columns[!is.na(columns)]
concord_cols <- concord_cols[grep(pattern = "missing", x = concord_cols, invert = T)]
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
concord_table_summary.df$prop_corr <- concord_table_summary.df$concord_count / (concord_table_summary.df$concord_count+ concord_table_summary.df$discord_count)

pdf(file = "05_compare/panel_vs_wgrs/hist_geno_concord_by_sample_prop_corr.pdf", width = 7.6, height = 4)
hist(x = concord_table_summary.df$prop_corr, las = 1, main = ""
     , xlab = "Proportion concordant"
     , breaks = 20
)
dev.off()

# Summary
summary(concord_table_summary.df$prop_corr)
sd(concord_table_summary.df$prop_corr, na.rm = T)


pdf(file = "05_compare/panel_vs_wgrs/ppn_concord_by_missing_data.pdf", width = 7.6, height = 4)
plot(x = concord_table_summary.df$missing_count, y = concord_table_summary.df$prop_corr
     , xlab = "Per sample missing data", ylab = "Proportion concordant"
     , las = 1
     )
dev.off()
