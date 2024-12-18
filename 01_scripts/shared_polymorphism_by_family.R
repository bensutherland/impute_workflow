# Assess shared polymorphism within each family, per family & per indiv
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
parent_vcf.FN    <- "06_screen_loci/parents_only_rename.vcf.gz" # parents, empirical
offspring_vcf.FN <- "06_screen_loci/offspring_only_rename.vcf.gz" # offspring, empirical

# Read in data
parent_vcf    <- read.vcfR(file = parent_vcf.FN)
parent_vcf
offspring_vcf <- read.vcfR(file = offspring_vcf.FN)
offspring_vcf


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


# Obtain vectors of loci names that are fully typed (no missing) per family
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
  labs(x = "Alt. allele freq.", y = NULL) +
  xlim(0,1)

F115_offspr_poly_parents_mono.plot <- F115_offspring_poly_parents_mono.af %>%
  ggplot(aes(x = alf2)) + 
  geom_histogram( bins = 40, fill = "darkgrey", color = "#e9ecef") +
  theme_bw() + 
  labs(x = "Alt. allele freq.", y = NULL) +
  xlim(0,1)

F115_offspr_poly_parents_poly.plot <- F115_offspring_poly_parents_poly.af %>%
  ggplot(aes(x = alf2)) + 
  geom_histogram( bins = 40, fill = "darkgrey", color = "#e9ecef") +
  theme_bw() + 
  labs(x = "Alt. allele freq.", y = NULL) +
  xlim(0,1)

F116_offspr_poly_parents_mono.plot <- F116_offspring_poly_parents_mono.af %>%
  ggplot(aes(x = alf2)) + 
  geom_histogram( bins = 40, fill = "darkgrey", color = "#e9ecef") +
  theme_bw() + 
  labs(x = "Alt. allele freq.", y = NULL) +
  xlim(0,1)

F116_offspr_poly_parents_poly.plot <- F116_offspring_poly_parents_poly.af %>%
  ggplot(aes(x = alf2)) + 
  geom_histogram( bins = 40, fill = "darkgrey", color = "#e9ecef") +
  theme_bw() + 
  labs(x = "Alt. allele freq.", y = NULL) +
  xlim(0,1)

F117_offspr_poly_parents_mono.plot <- F117_offspring_poly_parents_mono.af %>%
  ggplot(aes(x = alf2)) + 
  geom_histogram( bins = 40, fill = "darkgrey", color = "#e9ecef") +
  theme_bw() + 
  labs(x = "Alt. allele freq.", y = NULL) +
  xlim(0,1)

F117_offspr_poly_parents_poly.plot <- F117_offspring_poly_parents_poly.af %>%
  ggplot(aes(x = alf2)) + 
  geom_histogram( bins = 40, fill = "darkgrey", color = "#e9ecef") +
  theme_bw() + 
  labs(x = "Alt. allele freq.", y = NULL) +
  xlim(0,1)

# Multipanel plot
final.figure <- ggarrange(
  F114_offspr_poly_parents_poly.plot, F115_offspr_poly_parents_poly.plot, F116_offspr_poly_parents_poly.plot, F117_offspr_poly_parents_poly.plot
  , F114_offspr_poly_parents_mono.plot, F115_offspr_poly_parents_mono.plot, F116_offspr_poly_parents_mono.plot, F117_offspr_poly_parents_mono.plot
  , labels = c("A", "B", "C", "D", "E", "F", "G", "H")
  , ncol = 4, nrow = 2
  )

final.figure

pdf(file = "06_screen_loci/multipanel_AF_for_poly_off_and_par_and_poly_off_mono_par.pdf"
    , width = 9.6, height = 5
    )
print(final.figure)
dev.off()

# Save out results
# Remove unecessary items
rm(all_offspring.gl)
rm(parent_vcf)
rm(offspring_vcf)
rm(parent.gi)
rm(parent_gt.df)
gc()
#save.image(file = "tally_genos_per_sample/shared_polymorphs_up_to_indiv_offspring_analysis.RData")

# End


### Below is incomplete development, would need additional work ###

# #### 05. Offspring analysis, per individual ####
# # Make a list of the different offspring datasets
# # note: these have already been subset to only include loci that have *no missing data* in parents
# # Offspring data
# offspring.list <- list()
# offspring.list[["F114_offspring"]] <- F114_offspring.gl
# offspring.list[["F115_offspring"]] <- F115_offspring.gl
# offspring.list[["F116_offspring"]] <- F116_offspring.gl
# offspring.list[["F117_offspring"]] <- F117_offspring.gl
# 
# # Parent data
# parents.list <- list()
# parents.list[["F114_parents_nomissing"]] <- F114_parents_nomissing.gl
# parents.list[["F115_parents_nomissing"]] <- F115_parents_nomissing.gl
# parents.list[["F116_parents_nomissing"]] <- F116_parents_nomissing.gl
# parents.list[["F117_parents_nomissing"]] <- F117_parents_nomissing.gl
# 
# # Vectors of monomorphs in parental pairs
# parents_monomorphs.list <- list()
# parents_monomorphs.list[["F114_parents_monomorphs"]] <- F114_parents_monomorphs.vec
# parents_monomorphs.list[["F115_parents_monomorphs"]] <- F115_parents_monomorphs.vec
# parents_monomorphs.list[["F116_parents_monomorphs"]] <- F116_parents_monomorphs.vec
# parents_monomorphs.list[["F117_parents_monomorphs"]] <- F117_parents_monomorphs.vec
# 
# # Set names of fams
# families <- c("F114", "F115", "F116", "F117")
# 
# 
# # Set nulls
# foi <- NULL; indivs_for_family <- NULL
# offspring.oi.gl <- NULL; parent.oi.gl <- NULL; ind.oi.gl <- NULL; ind.oi_no_missing.gl <- NULL; ind.oi_no_missing_no_mono.gl <- NULL
# indiv <- NULL; indiv_polymorph_loci.vec <- NULL; parents_monomorphs.vec <- NULL; offspr_poly_parents_mono.vec <- NULL
# per_fam_info.df <- NULL
# per_fam_info.list <- list()
# parent.oi.gi <- NULL
# 
# # Loop over families
# for(i in 1:length(families)){
#   
#   # Select family of interest
#   foi <- families[i]
#   
#   # Select parent dataset for this family
#   parent.oi.gl <- parents.list[[grep(pattern = foi, x = names(parents.list))]]
#   print(paste0("Parents for ", foi, ":"))
#   print(indNames(parent.oi.gl))
#   parent.oi.gi <- gl2gi(parent.oi.gl) # convert to gi
#   pop(parent.oi.gi) <- rep(x = paste0(foi, "_parents"), times = nInd(parent.oi.gi))
#   
#   # Select offspring dataset for this family
#   offspring.oi.gl <- offspring.list[[grep(pattern = foi, x = names(offspring.list))]]
#   print("Offspring: ")
#   print(indNames(offspring.oi.gl))
#   offspring.oi.gi <- gl2gi(offspring.oi.gl)
#   pop(offspring.oi.gi) <- rep(x = paste0(foi, "_offspring"), times = nInd(offspring.oi.gi))
#   
#   # Empty the per_fam_info df and prepare it
#   per_fam_info.df <- data.frame(indiv=character(), 
#                                 num.poly=numeric(), 
#                                 num.mono=numeric(),
#                                 num.missing=numeric(),
#                                 num.poly.offspr.mono.parents=numeric(), 
#                                 stringsAsFactors=FALSE
#                                 ) 
#   
#   
#   for(o in 1:nInd(offspring.oi.gi)){
#     
#     # Subset to one offspring
#     ind.oi.gi <- offspring.oi.gi[o, ]
#     
#     # Retain the indiv name
#     per_fam_info.df[o,"indiv"] <- indNames(ind.oi.gi)
#     
#     # Combine the trio
#     trio.gl <- repool(ind.oi.gi, parent.oi.gi)
#     
#     pa <- private_alleles(gid = trio.gl)
#     
#     # Remove any missing loci from the indiv gl
#     ind.oi_no_missing.gl <- gl.filter.allna(ind.oi.gl)
#     
#     # Combine with parent data
#     
#     
#     
#     # Convert to genind
#     ind.oi.gi <- gl2gi(x = ind.oi.gl)
#     pop(in.oi.gi) <- "offspring"
#     
#     # Combine 
#     # combine above with parents, set pop, genind, determine private alleles
#     
#     
#     
#     
#     # Remove monomorphic
#     ind.oi_no_missing_no_mono.gl <- gl.filter.monomorphs(x = ind.oi_no_missing.gl, verbose = NULL)
#     # FAILS WITH SINGLE IND
#     
#     # Identify the missing loci
#     indiv_missing_loci.vec <- setdiff(x = locNames(ind.oi.gl), y = locNames(ind.oi_no_missing.gl))
#     
#     # Identify remaining loci as polymorphic
#     indiv_polymorph_loci.vec <- locNames(ind.oi_no_missing_no_mono.gl)
#     
#     # Identify the removed loci as monomorphic 
#     indiv_monomorph_loci.vec <- setdiff(x = locNames(ind.oi_no_missing.gl), y = locNames(ind.oi_no_missing_no_mono.gl))
#     
#     # Record how many polymorphic, monomorphic, and missing loci total for this indiv
#     per_fam_info.df[i,"num.poly"]    <- length(indiv_polymorph_loci.vec)
#     per_fam_info.df[i,"num.mono"]    <- length(indiv_monomorph_loci.vec)
#     per_fam_info.df[i,"num.missing"] <- length(indiv_missing_loci.vec)
#     
#     # How many of the polymorphic loci were monomorphic in parents? 
#     parents_monomorphs.vec <- parents_monomorphs.list[[grep(pattern = foi, x = names(parents_monomorphs.list))]]
#     offspr_poly_parents_mono.vec <- setdiff(x = indiv_polymorph_loci.vec, y = parents_monomorphs.vec)
#     per_fam_info.df[i,"num.poly.offspr.mono.parents"] <- length(offspr_poly_parents_mono.vec)
#     
#     # Retain info in list 
#     per_fam_info.list[[foi]] <- per_fam_info.df
#     
#   }
#   
# }
# 
# # Note: the above will miss alternate homozygotes, which could be captured by parent x offspring gid and private_alleles of poppr
