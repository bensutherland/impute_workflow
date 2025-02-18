# Summarize proportions or counts of each genotype per individual
# Adapted from Johnson et al. 2024 (G3) Population genomics of northern pike for impute_workflow
# B. Sutherland (2024-10-29)

# Clear space
#rm(list=ls())

# Load packages
#install.packages("vcfR")
#install.packages("rstudioapi")
#install.packages("dplyr")
#install.packages("tidyr")
#install.packages("ggplot2")
#install.packages("ggpubr")
library("tidyr")
library("vcfR")
library("rstudioapi")
library("dplyr")
library("ggplot2")
library("ggpubr")

# Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path)
setwd(current.path)
rm(current.path)

## Info
# sessionInfo()
options(scipen = 99999999)

#### 01. Set up ####
# Set filenames
vcf_empirical.FN   <- "06_screen_loci/parents_and_offspring_wgrs.vcf.gz" # all inds, empirical
vcf_imputed_ai2.FN <- "06_screen_loci/all_inds_wgrs_and_panel_biallele_only_ai2_imputed.vcf.gz" # all inds, imputed ai2
vcf_imputed_fi3.FN <- "06_screen_loci/all_inds_wgrs_and_panel_biallele_only_fi3_imputed.vcf.gz" # all inds, imputed fi3

# Select from the above three items
vcf.FN <- vcf_empirical.FN

# Read in data
vcf <- read.vcfR(file = vcf.FN)
vcf

# Extract genotypes
gt.df <- extract.gt(x = vcf, element = "GT") # note: do not use the 'as.numeric' method because there are more than two alleles in the VCF and this causes issues
gt.df[1:5,1:5]

# View all the sample names
colnames(gt.df)
colnames(gt.df) <- gsub(pattern = "ASY2-", replacement = "", x = colnames((gt.df)))
dim(gt.df)

# Manual inspection, test
table(gt.df[,1], useNA = "ifany")

# Inspect
gt.df <- as.data.frame(gt.df)

## Proof of principle
# test.df <- as.data.frame(table(gt.df$`Eluc-CL-M`, useNA = "ifany"))
# test.df

# Create summary df that collects the output of table in a long-form with multiple rows per sample, for each geno
result <- NULL; summary.df <- NULL; summary_all.df <- NULL
for(i in 1:ncol(gt.df)){
  
  # Debugging
  print(i)
  
  # Summarize using table
  result <- as.data.frame(table(gt.df[,i], useNA = "ifany"))
  
  # Create identifier
  name.vec <- rep(colnames(gt.df)[i], times = nrow(result))
  
  # Combine identifier and data
  summary_indiv.df <- cbind(name.vec, result)
  
  # Combine with other samples
  summary_all.df   <- rbind(summary_indiv.df, summary_all.df)
  
}

summary_all.df[1:15,]
summary_all.df <- as.data.frame(summary_all.df)
summary_all.df[1:15,]

colnames(summary_all.df)
colnames(summary_all.df) <- c("indiv", "geno", "count")
head(summary_all.df)

# # Restrict to only the target genos for plotting
# summary_all.df <- summary_all.df[which(summary_all.df$geno=="0/1" | 
#                                          summary_all.df$geno=="1/1"), ]

### Plot ###
# Reorder attempt
p <- ggplot(data = summary_all.df, aes(fill=geno, x = indiv, y = count)) +
  geom_bar(position='stack', stat = 'identity') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
        #, panel.border = element_blank()
        , panel.grid.major = element_blank()
        , panel.grid.minor = element_blank()
        #,
  ) +
  scale_fill_grey(start = 0.3, end = 0.7) + # grey colour
  xlab("Individual") + 
  ylab("Number variants") 
p

# Save plot
output_base.FN <- gsub(pattern = ".*\\/", replacement = "", vcf.FN)
output_base.FN <- gsub(pattern = "\\.vcf\\.gz", replacement = "", output_base.FN)
pdf(file = paste0("06_screen_loci/", output_base.FN, "_barplot_genos.pdf"), width = 9, height = 3.5)
print(p)
dev.off()

variants_per_indiv.plot <- p

# Calculate per population summary stats
head(summary_all.df)

sumstats.df  <- as.data.frame(summary_all.df)
head(sumstats.df)

#aggregate(x = sumstats.df$count, by = list(target = sumstats.df$indiv, sumstats.df$geno), FUN = sum)


# Separate the genos into multiple df then combine
homo_ref.df   <- sumstats.df[sumstats.df$geno %in% "0/0", ]
colnames(homo_ref.df)[colnames(homo_ref.df)=="count"] <- "count.homo_ref"
het.df        <- sumstats.df[sumstats.df$geno %in% "0/1", ]
colnames(het.df)[colnames(het.df)=="count"] <- "count.het"
homo_alt.df   <- sumstats.df[sumstats.df$geno %in% "1/1", ]
colnames(homo_alt.df)[colnames(homo_alt.df)=="count"] <- "count.homo_alt"
missing.df    <- sumstats.df[is.na(sumstats.df$geno), ]
colnames(missing.df)[colnames(missing.df)=="count"] <- "count.missing"

# Combine back together
wide.df <- merge(x = homo_ref.df, y = het.df, by = "indiv", all = T)
wide.df <- merge(x = wide.df, y = homo_alt.df, by = "indiv", all = T)
wide.df <- merge(x = wide.df, y = missing.df, by = "indiv", all = T)
# expect warnings()
head(wide.df)

# Select only needed
wide.df <- wide.df[,c("indiv", "count.homo_ref", "count.het", "count.homo_alt", "count.missing")]

# Checking
rowSums(wide.df[,2:ncol(wide.df)], na.rm = T)

write.table(x = wide.df, file = paste0("06_screen_loci/", output_base.FN, "_tallies.txt"), quote = F, sep = "\t", row.names = F)

# Save per sample genotype count plot as object
save(variants_per_indiv.plot, file = paste0("06_screen_loci/", output_base.FN, "_output.Rdata"))


#### DOWNSTREAM MULTI-PANEL PLOT ####
# To run the following multipanel plotting processes, you first need to run the above and generate Rdata for
# (1) all inds, empirical; (3) all inds, imputed by ai2; (4) all inds, imputed by fi3 

list.files(path = "06_screen_loci/", pattern = ".Rdata")

# Empirical data
load(file = "06_screen_loci/parents_and_offspring_wgrs_output.Rdata")
empirical.plot <- variants_per_indiv.plot
empirical.plot <- empirical.plot + theme(legend.position="none")

# Imputed data
load(file = "06_screen_loci/all_inds_wgrs_and_panel_biallele_only_ai2_imputed_output.Rdata")
all_inds_imputed_ai2.plot <- variants_per_indiv.plot
all_inds_imputed_ai2.plot <- all_inds_imputed_ai2.plot + theme(legend.position="none")

load(file = "06_screen_loci/all_inds_wgrs_and_panel_biallele_only_fi3_imputed_output.Rdata")
all_inds_imputed_fi3.plot <- variants_per_indiv.plot
all_inds_imputed_fi3.plot <- all_inds_imputed_fi3.plot + theme(legend.position="none")

multipanel.plot <- ggarrange(empirical.plot
                                          , all_inds_imputed_ai2.plot
                                          , all_inds_imputed_fi3.plot
                                          , ncol = 1, nrow = 3
                                          , labels = c("A", "B", "C")
)
multipanel.plot


pdf(file = "06_screen_loci/multipanel_geno_tally_per_sample.pdf", width = 20, height = 6)
print(multipanel.plot)
dev.off()

# Complete 
