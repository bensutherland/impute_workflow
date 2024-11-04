# Assess shared polymorphism within each family
# B. Sutherland (2024-11-01)

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
library("adegenet")
library("dartR")

# Set working directory to the ms_scallop_popgen repo
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path)
setwd(current.path)
rm(current.path)

## Info
# sessionInfo()
options(scipen = 99999999)

#### 01. Set up ####
# Set filenames
parent_vcf.FN    <- "tally_genos_per_sample/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_parents_only_rename.vcf.gz" # parents, empirical
offspring_vcf.FN <- "tally_genos_per_sample/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_offspring_only_rename.vcf.gz" # offspring, empirical

# Read in data
parent_vcf    <- read.vcfR(file = parent_vcf.FN)
offspring_vcf <- read.vcfR(file = offspring_vcf.FN)

#### 02. Determine which loci are fully typed in each parental pairing ####
# Extract genotypes
parent_gt.df <- extract.gt(x = parent_vcf, element = "GT")
dim(parent_gt.df)
parent_gt.df[1:5,1:5]
colnames(parent_gt.df)

# Subset to families
# F114 parents
F114_parents.df <- parent_gt.df[,c("65-4F", "58-33M")] # select the parents for this family
colSums(is.na(F114_parents.df))                        # how many NAs per parent?
F114_parents.df <- as.data.frame(F114_parents.df)      # make df
F114_parents.df$num.missing <- rowSums(is.na(F114_parents.df)) # count NAs per locus
F114_parents.df <- F114_parents.df[F114_parents.df$num.missing == 0, ] # keep only those loci with no missing data
dim(F114_parents.df)  # how many remain? 

# F115 parents
F115_parents.df <- parent_gt.df[,c("65-8F", "79-13M")] # select the parents for this family
colSums(is.na(F115_parents.df))                        # how many NAs per parent?
F115_parents.df <- as.data.frame(F115_parents.df)      # make df
F115_parents.df$num.missing <- rowSums(is.na(F115_parents.df)) # count NAs per locus
F115_parents.df <- F115_parents.df[F115_parents.df$num.missing == 0, ] # keep only those loci with no missing data
dim(F115_parents.df)  # how many remain? 

# F116 parents
F116_parents.df <- parent_gt.df[,c("55-41F", "65-19M")] # select the parents for this family
colSums(is.na(F116_parents.df))                        # how many NAs per parent?
F116_parents.df <- as.data.frame(F116_parents.df)      # make df
F116_parents.df$num.missing <- rowSums(is.na(F116_parents.df)) # count NAs per locus
F116_parents.df <- F116_parents.df[F116_parents.df$num.missing == 0, ] # keep only those loci with no missing data
dim(F116_parents.df)  # how many remain? 

# F117 parents
F117_parents.df <- parent_gt.df[,c("58-9F", "79-1M")] # select the parents for this family
colSums(is.na(F117_parents.df))                        # how many NAs per parent?
F117_parents.df <- as.data.frame(F117_parents.df)      # make df
F117_parents.df$num.missing <- rowSums(is.na(F117_parents.df)) # count NAs per locus
F117_parents.df <- F117_parents.df[F117_parents.df$num.missing == 0, ] # keep only those loci with no missing data
dim(F117_parents.df)  # how many remain? 


# Identify per family which loci are fully typed (no NAs)
F114_parents_loci_no_missing.vec <- rownames(F114_parents.df)
F114_parents_loci_no_missing.vec <- gsub(pattern = "\\.", replacement = "_", x = F114_parents_loci_no_missing.vec)
length(F114_parents_loci_no_missing.vec)
F115_parents_loci_no_missing.vec <- rownames(F115_parents.df)
F115_parents_loci_no_missing.vec <- gsub(pattern = "\\.", replacement = "_", x = F115_parents_loci_no_missing.vec)
length(F115_parents_loci_no_missing.vec)
F116_parents_loci_no_missing.vec <- rownames(F116_parents.df)
F116_parents_loci_no_missing.vec <- gsub(pattern = "\\.", replacement = "_", x = F116_parents_loci_no_missing.vec)
length(F116_parents_loci_no_missing.vec)
F117_parents_loci_no_missing.vec <- rownames(F117_parents.df)
F117_parents_loci_no_missing.vec <- gsub(pattern = "\\.", replacement = "_", x = F117_parents_loci_no_missing.vec)
length(F117_parents_loci_no_missing.vec)


#### 03. Subset parental data to only fully typed loci ####
parent.gi <- vcfR2genind(x = parent_vcf, sep = "/", return.alleles = F)

# Subset to specific parents
# F114 
F114_parents.gi <- parent.gi[i = c("65-4F", "58-33M")]
F114_parents.gl <- gi2gl(F114_parents.gi)
rm(F114_parents.gi)

# F115 
F115_parents.gi <- parent.gi[i = c("65-8F", "79-13M")]
F115_parents.gl <- gi2gl(F115_parents.gi)
rm(F115_parents.gi)

# F116 
F116_parents.gi <- parent.gi[i = c("55-41F", "65-19M")]
F116_parents.gl <- gi2gl(F116_parents.gi)
rm(F116_parents.gi)

# F117 
F117_parents.gi <- parent.gi[i = c("58-9F", "79-1M")]
F117_parents.gl <- gi2gl(F117_parents.gi)
rm(F117_parents.gi)

# Subset to only keep the no-missing loci
# F114
F114_parents_nomissing.gl <- F114_parents.gl[, (locNames(F114_parents.gl) %in% F114_parents_loci_no_missing.vec)]
F114_parents_nomissing.gl

# F115
F115_parents_nomissing.gl <- F115_parents.gl[, (locNames(F115_parents.gl) %in% F115_parents_loci_no_missing.vec)]
F115_parents_nomissing.gl

# F116
F116_parents_nomissing.gl <- F116_parents.gl[, (locNames(F116_parents.gl) %in% F116_parents_loci_no_missing.vec)]
F116_parents_nomissing.gl

# F117
F117_parents_nomissing.gl <- F117_parents.gl[, (locNames(F117_parents.gl) %in% F117_parents_loci_no_missing.vec)]
F117_parents_nomissing.gl

#### 04. Identify polymorphic and monomorphic loci in parental pairs ####
# Find polymorphic and monomorphic loci in parents, make separate lists. Per family. 
# F114
F114_parents_nomissing_nomonomorphs.gl <- gl.filter.monomorphs(x = F114_parents_nomissing.gl, verbose = NULL)
F114_parents_polymorphs.vec <- locNames(F114_parents_nomissing_nomonomorphs.gl)
length(F114_parents_polymorphs.vec)
F114_parents_monomorphs.vec <- setdiff(x = locNames(F114_parents_nomissing.gl), y = F114_parents_polymorphs.vec)
length(F114_parents_monomorphs.vec)

# F115
F115_parents_nomissing_nomonomorphs.gl <- gl.filter.monomorphs(x = F115_parents_nomissing.gl, verbose = NULL)
F115_parents_polymorphs.vec <- locNames(F115_parents_nomissing_nomonomorphs.gl)
length(F115_parents_polymorphs.vec)
F115_parents_monomorphs.vec <- setdiff(x = locNames(F115_parents_nomissing.gl), y = F115_parents_polymorphs.vec)
length(F115_parents_monomorphs.vec)

# F116
F116_parents_nomissing_nomonomorphs.gl <- gl.filter.monomorphs(x = F116_parents_nomissing.gl, verbose = NULL)
F116_parents_polymorphs.vec <- locNames(F116_parents_nomissing_nomonomorphs.gl)
length(F116_parents_polymorphs.vec)
F116_parents_monomorphs.vec <- setdiff(x = locNames(F116_parents_nomissing.gl), y = F116_parents_polymorphs.vec)
length(F116_parents_monomorphs.vec)

# F117
F117_parents_nomissing_nomonomorphs.gl <- gl.filter.monomorphs(x = F117_parents_nomissing.gl, verbose = NULL)
F117_parents_polymorphs.vec <- locNames(F117_parents_nomissing_nomonomorphs.gl)
length(F117_parents_polymorphs.vec)
F117_parents_monomorphs.vec <- setdiff(x = locNames(F117_parents_nomissing.gl), y = F117_parents_polymorphs.vec)
length(F117_parents_monomorphs.vec)


#### 04. Offspring analysis ####
# Convert from vcfR to genlight
all_offspring.gl <- vcfR2genlight(offspring_vcf) # convert to genlight
locNames(all_offspring.gl) <- gsub(pattern = "\\.", replacement = "_", x = locNames(all_offspring.gl)) # format to match

# Subset by family
F114_offspring.gl <- all_offspring.gl[grep(pattern = "-114-", x = indNames(all_offspring.gl)), ]
F114_offspring.gl

F115_offspring.gl <- all_offspring.gl[grep(pattern = "-115-", x = indNames(all_offspring.gl)), ]
F115_offspring.gl

F116_offspring.gl <- all_offspring.gl[grep(pattern = "-116-", x = indNames(all_offspring.gl)), ]
F116_offspring.gl

F117_offspring.gl <- all_offspring.gl[grep(pattern = "-117-", x = indNames(all_offspring.gl)), ]
F117_offspring.gl

# Only keep loci that are 'no.missing' from their respective parents
F114_offspring.gl <- F114_offspring.gl[, (locNames(F114_offspring.gl) %in% F114_parents_loci_no_missing.vec)]
F114_offspring.gl

F115_offspring.gl <- F115_offspring.gl[, (locNames(F115_offspring.gl) %in% F115_parents_loci_no_missing.vec)]
F115_offspring.gl

F116_offspring.gl <- F116_offspring.gl[, (locNames(F116_offspring.gl) %in% F116_parents_loci_no_missing.vec)]
F116_offspring.gl

F117_offspring.gl <- F117_offspring.gl[, (locNames(F117_offspring.gl) %in% F117_parents_loci_no_missing.vec)]
F117_offspring.gl


# # OPTIONAL # Remove missing data (>70% per family) [ probably not needed ]
# F114_offspring.gl <- gl.filter.callrate(x = F114_offspring.gl, method = "loc", threshold = 0.70, mono.rm = F)
# F114_offspring.gl

# Identify polymorphs and monomorphs in offspring
F114_offspring_polymorphs.gl <- gl.filter.monomorphs(x = F114_offspring.gl, verbose = NULL)
F114_offspring_polymorphs.vec <- locNames(F114_offspring_polymorphs.gl)
length(F114_offspring_polymorphs.vec)
F114_offspring_monomorphs.vec <- setdiff(x = locNames(F114_offspring.gl), y = F114_offspring_polymorphs.vec)
length(F114_offspring_monomorphs.vec)

F115_offspring_polymorphs.gl <- gl.filter.monomorphs(x = F115_offspring.gl, verbose = NULL)
F115_offspring_polymorphs.vec <- locNames(F115_offspring_polymorphs.gl)
length(F115_offspring_polymorphs.vec)
F115_offspring_monomorphs.vec <- setdiff(x = locNames(F115_offspring.gl), y = F115_offspring_polymorphs.vec)
length(F115_offspring_monomorphs.vec)

F116_offspring_polymorphs.gl <- gl.filter.monomorphs(x = F116_offspring.gl, verbose = NULL)
F116_offspring_polymorphs.vec <- locNames(F116_offspring_polymorphs.gl)
length(F116_offspring_polymorphs.vec)
F116_offspring_monomorphs.vec <- setdiff(x = locNames(F116_offspring.gl), y = F116_offspring_polymorphs.vec)
length(F116_offspring_monomorphs.vec)

F117_offspring_polymorphs.gl <- gl.filter.monomorphs(x = F117_offspring.gl, verbose = NULL)
F117_offspring_polymorphs.vec <- locNames(F117_offspring_polymorphs.gl)
length(F117_offspring_polymorphs.vec)
F117_offspring_monomorphs.vec <- setdiff(x = locNames(F117_offspring.gl), y = F117_offspring_polymorphs.vec)
length(F117_offspring_monomorphs.vec)

# How many loci are polymorphic in offspring and monomorphic in parents? 
F114_offspring_poly_parents_mono.vec <- intersect(x = F114_offspring_polymorphs.vec, y = F114_parents_monomorphs.vec)
length(F114_offspring_poly_parents_mono.vec)
length(F114_offspring_poly_parents_mono.vec) / length(F114_offspring_polymorphs.vec)

F114_offspring_poly_parents_poly.vec <- intersect(x = F114_offspring_polymorphs.vec, y = F114_parents_polymorphs.vec)
length(F114_offspring_poly_parents_poly.vec)
length(F114_offspring_poly_parents_poly.vec) / length(F114_offspring_polymorphs.vec)

F115_offspring_poly_parents_mono.vec <- intersect(x = F115_offspring_polymorphs.vec, y = F115_parents_monomorphs.vec)
length(F115_offspring_poly_parents_mono.vec)
length(F115_offspring_poly_parents_mono.vec) / length(F115_offspring_polymorphs.vec)

F115_offspring_poly_parents_poly.vec <- intersect(x = F115_offspring_polymorphs.vec, y = F115_parents_polymorphs.vec)
length(F115_offspring_poly_parents_poly.vec)
length(F115_offspring_poly_parents_poly.vec) / length(F115_offspring_polymorphs.vec)

F116_offspring_poly_parents_mono.vec <- intersect(x = F116_offspring_polymorphs.vec, y = F116_parents_monomorphs.vec)
length(F116_offspring_poly_parents_mono.vec)
length(F116_offspring_poly_parents_mono.vec) / length(F116_offspring_polymorphs.vec)

F116_offspring_poly_parents_poly.vec <- intersect(x = F116_offspring_polymorphs.vec, y = F116_parents_polymorphs.vec)
length(F116_offspring_poly_parents_poly.vec)
length(F116_offspring_poly_parents_poly.vec) / length(F116_offspring_polymorphs.vec)

F117_offspring_poly_parents_mono.vec <- intersect(x = F117_offspring_polymorphs.vec, y = F117_parents_monomorphs.vec)
length(F117_offspring_poly_parents_mono.vec)
length(F117_offspring_poly_parents_mono.vec) / length(F117_offspring_polymorphs.vec)

F117_offspring_poly_parents_poly.vec <- intersect(x = F117_offspring_polymorphs.vec, y = F117_parents_polymorphs.vec)
length(F117_offspring_poly_parents_poly.vec)
length(F117_offspring_poly_parents_poly.vec) / length(F117_offspring_polymorphs.vec)


# Calculate AF of loci that are poly in offspring and mono in parents, or poly in both 
F114_offspring_poly_parents_mono.af <- gl.alf(x = F114_offspring.gl[, F114_offspring_poly_parents_mono.vec]) # poly are mono in parents
F114_offspring_poly_parents_poly.af <- gl.alf(x = F114_offspring.gl[, F114_offspring_poly_parents_poly.vec]) # all polymorphs

F115_offspring_poly_parents_mono.af <- gl.alf(x = F115_offspring.gl[, F115_offspring_poly_parents_mono.vec]) # poly are mono in parents
F115_offspring_poly_parents_poly.af <- gl.alf(x = F115_offspring.gl[, F115_offspring_poly_parents_poly.vec]) # all polymorphs

F116_offspring_poly_parents_mono.af <- gl.alf(x = F116_offspring.gl[, F116_offspring_poly_parents_mono.vec]) # poly are mono in parents
F116_offspring_poly_parents_poly.af <- gl.alf(x = F116_offspring.gl[, F116_offspring_poly_parents_poly.vec]) # all polymorphs

F117_offspring_poly_parents_mono.af <- gl.alf(x = F117_offspring.gl[, F117_offspring_poly_parents_mono.vec]) # poly are mono in parents
F117_offspring_poly_parents_poly.af <- gl.alf(x = F117_offspring.gl[, F117_offspring_poly_parents_poly.vec]) # all polymorphs


# Plot histogram
# par(mfrow = c(2,4), mar = c(5.1, 4.1, 4.1, 2.1))
# hist(F114_offspring_poly_parents_poly.af[,"alf2"], main = "F114, offspr. poly, parents poly"
#      , xlab = "Allele Freq. (alt allele)"
#      , ylab = ""
#      , las = 1) # Plot AF of alt allele
# hist(F114_offspring_poly_parents_mono.af[,"alf2"], main = "F114, offspr. poly, parents mono"
#      , xlab = "Allele Freq. (alt allele)"
#      , ylab = ""
#      , las = 1) # Plot AF of alt allele

F114_offspr_poly_parents_mono.plot <- F114_offspring_poly_parents_mono.af %>%
  ggplot(aes(x = alf2)) + 
  geom_histogram( bins = 40, fill = "darkgrey", color = "#e9ecef") +
  theme_bw() + 
  labs(x = "Alt. allele freq.", y = NULL)
# F114_offspr_poly_parents_mono.plot <- F114_offspr_poly_parents_mono.plot + annotate("text", x = 0.65, y = 250000
#                                               , label = paste0(length(F114_offspring_poly_parents_mono.vec), " loci")
#                                       )


F114_offspr_poly_parents_poly.plot <- F114_offspring_poly_parents_poly.af %>%
  ggplot(aes(x = alf2)) + 
  geom_histogram( bins = 40, fill = "darkgrey", color = "#e9ecef") +
  theme_bw() + 
  labs(x = "Alt. allele freq.", y = NULL)

F115_offspr_poly_parents_mono.plot <- F115_offspring_poly_parents_mono.af %>%
  ggplot(aes(x = alf2)) + 
  geom_histogram( bins = 40, fill = "darkgrey", color = "#e9ecef") +
  theme_bw() + 
  labs(x = "Alt. allele freq.", y = NULL)

F115_offspr_poly_parents_poly.plot <- F115_offspring_poly_parents_poly.af %>%
  ggplot(aes(x = alf2)) + 
  geom_histogram( bins = 40, fill = "darkgrey", color = "#e9ecef") +
  theme_bw() + 
  labs(x = "Alt. allele freq.", y = NULL)

F116_offspr_poly_parents_mono.plot <- F116_offspring_poly_parents_mono.af %>%
  ggplot(aes(x = alf2)) + 
  geom_histogram( bins = 40, fill = "darkgrey", color = "#e9ecef") +
  theme_bw() + 
  labs(x = "Alt. allele freq.", y = NULL)

F116_offspr_poly_parents_poly.plot <- F116_offspring_poly_parents_poly.af %>%
  ggplot(aes(x = alf2)) + 
  geom_histogram( bins = 40, fill = "darkgrey", color = "#e9ecef") +
  theme_bw() + 
  labs(x = "Alt. allele freq.", y = NULL)

F117_offspr_poly_parents_mono.plot <- F117_offspring_poly_parents_mono.af %>%
  ggplot(aes(x = alf2)) + 
  geom_histogram( bins = 40, fill = "darkgrey", color = "#e9ecef") +
  theme_bw() + 
  labs(x = "Alt. allele freq.", y = NULL)

F117_offspr_poly_parents_poly.plot <- F117_offspring_poly_parents_poly.af %>%
  ggplot(aes(x = alf2)) + 
  geom_histogram( bins = 40, fill = "darkgrey", color = "#e9ecef") +
  theme_bw() + 
  labs(x = "Alt. allele freq.", y = NULL)

# Multipanel plot
final.figure <- ggarrange(
  F114_offspr_poly_parents_poly.plot, F115_offspr_poly_parents_poly.plot, F116_offspr_poly_parents_poly.plot, F117_offspr_poly_parents_poly.plot
  , F114_offspr_poly_parents_mono.plot, F115_offspr_poly_parents_mono.plot, F116_offspr_poly_parents_mono.plot, F117_offspr_poly_parents_mono.plot
  , labels = c("A", "B", "C", "D", "E", "F", "G", "H")
  , ncol = 4, nrow = 2
  )

final.figure

pdf(file = "tally_genos_per_sample/multipanel_AF_for_poly_off_and_par_and_poly_off_mono_par.pdf"
    , width = 9.6, height = 5
    )
print(final.figure)
dev.off()


# #### OLD CODE ####
# # Subset to individual
# test <- parent_vcf[j = "55-41F"]
# 
# 
# 
# 
# # Subset to individual families (parents)
# all_parent.gl <- vcfR2genlight(parent_vcf) # convert to genlight
# indNames(all_parent.gl) # identify inds
# pop(all_parent.gl) <- c("F116", "F114", "F116", "F114", "F115", "F117", "F115", "F117") # assign pop IDs
# separated_pops_parents.list <- seppop(all_parent.gl) # separate to indiv pops
# F114_parents.gl <- separated_pops_parents.list$F114
# F114_parents.gl$gen[[1]]@snp # gives per-sample missing data
# #F114_parents.gi <- gl2gi(x = F114_parents.gl) # throws error
# 
# 
# 
# 
# # Identify fully typed loci
# F114_parents.gl <- gl.filter.monomorphs(x = F114_parents.gl)
# 
# 
