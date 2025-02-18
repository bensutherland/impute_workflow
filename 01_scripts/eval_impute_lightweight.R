# Read in imputation outputs in ai2 format and compare with empirical data (ai2 format) of same individuals to evaluate imputation
# Importantly, keep in mind that if working with ai2, 9 or 5 = missing for ai2 and fi3, respectively
# Currently assumes that there is no missing data in imputed, and that any missing data in empirical is coded as 9

# B. Sutherland (2024-07-25)

### Front Matter ####
# Clean space
# rm(list=ls())

## Install and load packages
#install.packages("rstudioapi")
#install.packages("data.table")
library("rstudioapi")
library("data.table")
library(tidyr)
library(dplyr)
library("matrixStats")

# Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "/01_scripts", replacement = "", x = current.path)
setwd(current.path)
rm(current.path)

## Info
# sessionInfo()

options(scipen = 9999999)

# Set variables
offspring_imputed_ai2.FN     <- "05_compare/all_chr_combined.txt" # imputed ai2 FN
offspring_imputed_fi3.FN <- "04_impute/fimpute/fi3_loci_by_inds_all_imputed_chr.txt" # impute fi3 FN
offspring_10X_ai2.FN         <- "05_compare/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_offspring_only_rename_ai2.txt" # 10X genotypes
subsampled <- FALSE # TRUE or FALSE whether this is subsampled data from all high density
output_folder <- "05_compare"

# Select the imputation software target to use
impute_software <- "fi3" # fi3 or ai2


#### 01. Load data ####
# Read in imputed data
if(impute_software=="ai2"){
  
  # Reporting
  print(paste0("Reading in data imputed by ", impute_software))
  
  imputed.df <- fread(file = offspring_imputed_ai2.FN, sep = "\t")
  dim(imputed.df)
  imputed.df <- as.data.frame(imputed.df) # convert to df
  imputed.df[1:5,1:5]

}else if(impute_software=="fi3"){
  
  # Reporting
  print(paste0("Reading in data imputed by ", impute_software))
  
  imputed.df <- fread(file = offspring_imputed_fi3.FN, sep = "\t")
  dim(imputed.df)
  imputed.df <- as.data.frame(imputed.df)
  colnames(imputed.df)[grep(pattern = "V1", x = colnames(imputed.df))] <- "mname"
  imputed.df$mname <- gsub(pattern = "__", replacement = " ", x = imputed.df$mname) # convert to mname as per ai2
  imputed.df[1:5,1:5]

}

# Reporting
print(paste0("Number loci in imputed dataset: ", nrow(imputed.df)))
print(paste0("Number indivs in imputed dataset: ", length(grep(pattern = "mname", x = colnames(imputed.df), invert = T))))

# Read in empirical data (10X)
print(paste0("Reading in empirical data"))
empirical.df <- fread(file = offspring_10X_ai2.FN, sep = "\t")
dim(empirical.df)
empirical.df <- as.data.frame(empirical.df) # convert to df
empirical.df[1:5,1:5]

# Reporting
print(paste0("Number loci: ", nrow(empirical.df)))
print(paste0("Number indivs in empirical dataset: ", length(grep(pattern = "mname", x = colnames(empirical.df), invert = T))))

#### 01.2. Limit to offspring only ####
## Make sample names match for imputed data
head(colnames(imputed.df), n = 20)

# Only keep offspring
imputed.df <- imputed.df[, grep(pattern = "mname|ASY2", x = colnames(imputed.df))] # remove parents, keep mname and ASY2 inds
dim(imputed.df)


#### 01.3. If working with subset, remove the training from the test ####
# # Remove 10X data (if present; this was used to support imputation in one evaluation)
# imputed.df <- imputed.df[, grep(pattern = "fastq.gz", x = colnames(imputed.df), invert = T)] # remove parents, keep mname and ASY2 inds
# dim(imputed.df)


#### 02. Prepare data for matching ####
# Check for duplicates (stop if there are any)
table(duplicated(colnames(imputed.df))) # any duplicates?

# If working with subsampled data, rather than panel data
if(subsampled == TRUE){
  
  # Remove everything after the first underscore
  colnames.df <- gsub(pattern = "_.*", replacement = "", x = colnames(imputed.df))
  colnames.df <- as.data.frame(colnames.df)
  colnames.df <- colnames.df[2:nrow(colnames.df), ] # drop mname for now
  colnames.df <- as.data.frame(colnames.df)
  colnames.df <- separate(data = colnames.df, col = "colnames.df", into = c("assay", "fam", "rep", "ind", "noise"), sep = "-", remove = T)
  colnames.df$name <- paste0(colnames.df$assay, "_", colnames.df$fam, "_", colnames.df$rep, "_", colnames.df$ind)
  colnames(imputed.df) <- c("mname", colnames.df$name)
  
}else{
  
  print("Not working with subsampled data.")
  
}

## Make sample names match for empirical data
# /NOTE: no longer doing this, samples should already match #
# /NOTE: also removes requirement for checking for duplicate samples, these are dropped earlier #
head(colnames(empirical.df))
head(colnames(imputed.df))


#### 03. Matching ####
## What columns match between the two
keep.cols <- intersect(x = colnames(imputed.df), y = colnames(empirical.df))
length(keep.cols)

# Keep only the keep cols from imputed
imputed.df <- imputed.df[, colnames(imputed.df) %in% keep.cols]
dim(imputed.df)

if(subsampled==TRUE){
  
  drop_cols <- c("ASY2_117_R2_7.1", "ASY2_117_R5_9.1")
  
  imputed.df <- imputed.df[, !(colnames(imputed.df) %in% drop_cols) ]
  dim(imputed.df)
  
}

# Keep only the keep cols from empirical
empirical.df <- empirical.df[, colnames(empirical.df) %in% keep.cols]
dim(empirical.df)


#### 04. Clean up mname
imputed.df$mname <- gsub(pattern = " ", replacement = "__", x = imputed.df$mname)
imputed.df[1:5,1:5]
empirical.df$mname <- gsub(pattern = " ", replacement = "__", x = empirical.df$mname)
empirical.df[1:5,1:5]


#### 05. Add unique identifier to separate the two datasets
colnames(imputed.df)   <- paste0(colnames(imputed.df), "_imputed")
colnames(empirical.df) <- paste0(colnames(empirical.df), "_empirical")


#### 06. Merge the two datasets
all_data.df <- merge(x = imputed.df, y = empirical.df, by.x = "mname_imputed", by.y = "mname_empirical")
dim(all_data.df)
all_data.df[1:5,1:5]
colnames(all_data.df) <- gsub("mname_imputed", replacement = "mname", x = colnames(all_data.df)) # clean up mname colname

# Sort by colname, then put mname first
all_data.df <- all_data.df[ , order(colnames(all_data.df))]
all_data.df <- all_data.df %>% 
  select("mname", everything())
all_data.df[1:5,1:5]


# Create output foldername
date <- format(Sys.time(), "%Y-%m-%d")
output_folder <- paste0(output_folder, "/concord_eval_", impute_software, "_", date)
print(paste0("Making output directory: ", output_folder))
dir.create(path = output_folder, showWarnings = F)

#### SETUP CONCORD EVAL ####
# Identify chr
chr <- unique(gsub(pattern = "__.*", replacement = "", x = all_data.df$mname))
chr

# Identify samples
samples.vec <- unique(gsub(pattern = "_empirical|_imputed", replacement = "", x = colnames(all_data.df)))
samples.vec <- samples.vec[grep(pattern = "mname", x = samples.vec, invert = T)]


#### 04. Evaluate per-sample overall concordance ####
soi <- NULL ; score <- NULL

# Loop to evaluate concordance per sample
# Set up an empty df to fill
result.df <- matrix(data = NA, nrow = length(samples.vec), ncol = 5)
result.df <- as.data.frame(result.df)
colnames(result.df) <- c("sample", "num.loci", "num.match", "num.missing", "prop.match")

soi <- NULL ; score <- NULL; calc.cor <- NULL; imputed_vector <- NULL; empirical_vector <- NULL
pair_compare.df <- NULL

locus_info.df <- matrix(data = NA, nrow = nrow(all_data.df), ncol = length(samples.vec) + 1)
colnames(locus_info.df) <- c("mname", samples.vec)
locus_info.df <- as.data.frame(locus_info.df)
locus_info.df$mname <- all_data.df$mname
locus_info.df[1:5,1:5]
dim(locus_info.df)



# Loop for each individual
for(i in 1:length(samples.vec)){
  
  # Select the sample of interest
  soi <- samples.vec[i]
  
  # Select data for the sample of interest, both comparisons
  pair_compare.df <- all_data.df[, c(paste0(soi, "_imputed"), paste0(soi, "_empirical"))]
  head(pair_compare.df, n = 5)
  
  # Per locus comparison (match = TRUE); add result to the locus info file
  locus_info.df[, i+1] <- pair_compare.df[,1] == pair_compare.df[,2]
  #locus_info.df[1:5,1:5]
  
  # If the line is missing in pair_compare (1 or 2) make it an NA in the locus_info file
  locus_info.df[pair_compare.df[,1]=="9", i+1] <- NA # if there is a missing record, NA the line
  locus_info.df[pair_compare.df[,2]=="9", i+1] <- NA # if there is a missing record, NA the line
  
  ### TODO: confirm do not need to set as "5" for fi3
  
  # For summarizing, go back to pair_compare and drop any row that has missing data in either sample
  #  i.e., keep complete records only
  dim(pair_compare.df)
  pair_compare.df <- pair_compare.df[pair_compare.df[,1]!="9", ]
  pair_compare.df <- pair_compare.df[pair_compare.df[,2]!="9", ]
  
  # Per individual comparison
  # After have dropped any missing records, 
  #   sum up the number of identical matches between the two
  score <- sum(pair_compare.df[,1] == pair_compare.df[,2])
  
  # tally the number of missing values for this pair
  num_missing <- nrow(all_data.df) - nrow(pair_compare.df)
  
  # calculate the proportion correct for this sample by 
  # dividing the number correct by the total without any missing data
  prop_correct <- score / nrow(pair_compare.df)
  
  # Store results per sample
  result.df[i,"sample"]      <- soi
  result.df[i,"num.loci"]    <- nrow(all_data.df)
  result.df[i,"num.match"]   <- score
  result.df[i,"num.missing"] <- num_missing
  result.df[i,"prop.match"]  <- round(x = prop_correct, digits = 3)
  
}

head(result.df)
write.table(x = result.df, file = paste0(output_folder, "/per_sample_overall_concord.txt"), sep = "\t", row.names = F)

# Summarize
summary(result.df$num.match)
summary(result.df$prop.match)
sd(result.df$prop.match)

# Plots
# Plot proportion match by number missing in empirical
pdf(file = paste0(output_folder, "/per_sample_prop_match_by_number_missing.pdf"), width = 6, height = 5)
plot(x = result.df$num.missing, y = result.df$prop.match, las = 1)
dev.off()

# Plot proportion match by number missing in empirical
pdf(file = paste0(output_folder, "/per_sample_prop_match_histogram.pdf"), width = 6, height = 5)
hist(x = result.df$prop.match, main = "", xlab = "Proportion matching (empirical & imputed)", las = 1)
dev.off()


##### Assess per locus stats #####
### NEED TO REDO, create a new, separate df ####
# rownames(locus_info.df) <- locus_info.df$mname
# locus_info.df <- locus_info.df[,2:ncol(locus_info.df)]
# 
# # Sum the number of matches per locus
# locus_info.df$matches <- rowSums(x = locus_info.df, na.rm = T)
# # Sum the number of missing values per locus
# locus_info.df$NA_vals <- rowSums(is.na(locus_info.df))
# locus_info.df$num_comparisons <- length(samples.vec) - locus_info.df$NA_vals
# 
# locus_info.df$percent.corr <- NA
# 
# # # Identify the prop corr per marker
# # for(i in 1:nrow(locus_info.df)){
# #   
# #   locus_info.df$percent.corr[i] <- locus_info.df$matches[i] / (length(samples.vec) - locus_info.df$NA_vals[i])
# #   
# # }
# 
# 
# pdf(file = paste0(run_folder, "/concord_eval_per_locus_percent_corr_hist.pdf"), width = 6.5, height = 3.5)
# hist(x = locus_info.df$percent.corr, main = "", breaks = 10
#      , xlab = paste0("Per locus percent of complete records concordant")
#      , las = 1
# )
# #text(x = 0.67, y = 30, paste0(nrow(subset_data.df), " loci"))
# dev.off()



#### 05. Per chromosome, evaluate concordance ####

# Per chromosome
subset_data.df <- NULL; result.df <- NULL; soi <- NULL ; score <- NULL
for(c in 1:length(chr)){
  
  print(paste0("Working on chr ", chr[c]))
  
  # Subset to the target chr
  subset_data.df <- all_data.df[grep(pattern = chr[c], x = all_data.df$mname), ]
  print(paste0("Number of records in this chr: ", nrow(subset_data.df)))
  
  # Loop to evaluate concordance per sample
  # Set up an empty df to fill
  result.df <- matrix(data = NA, nrow = length(samples.vec), ncol = 6)
  result.df <- as.data.frame(result.df)
  colnames(result.df) <- c("sample", "num.loci", "num.match", "num.missing", "prop.match", "pearson.cor")
  
  soi <- NULL ; score <- NULL; calc.cor <- NULL; imputed_vector <- NULL; empirical_vector <- NULL
  for(i in 1:length(samples.vec)){
    
    # Select the sample of interest
    soi <- samples.vec[i]
    
    # NOTE: TODO: should convert "9" values to NAs here
    
    # sum up the number of identical matches between the empirical and the imputed for this sample
    score <- sum(subset_data.df[, paste0(soi, "_empirical")] == subset_data.df[, paste0(soi, "_imputed")])
    
    # tally the number of missing values in either of the subsets (note: the imputed will always be no missing)
    num_missing <- sum(subset_data.df[, paste0(soi, "_empirical")]==9) + sum(subset_data.df[, paste0(soi, "_imputed")]==9)
    #TODO: the logic of the above num_missing calc should be improved, as a simple sum of the two would not
    #  technically work if there were actually 9s (missing vals) in the imputed
    
    # calculate the proportion correct for this sample by dividing the number correct by the total (with the number missing subtracted)
    prop_correct <- score / (nrow(subset_data.df) - num_missing)
    
    # Calculate the pearson correlation between the two datatypes, only consider complete obs
    imputed_vector   <- subset_data.df[, paste0(soi, "_imputed")]
    empirical_vector <- subset_data.df[, paste0(soi, "_empirical")]
    imputed_vector   <- gsub(pattern = "9", replacement = NA, x = imputed_vector)
    empirical_vector <- gsub(pattern = "9", replacement = NA, x = empirical_vector)
    
    imputed_vector <- as.numeric(imputed_vector)
    empirical_vector <- as.numeric(empirical_vector)
    
    calc.cor <-  cor(x = imputed_vector, y = empirical_vector
                   , method = "pearson", use = "pairwise.complete.obs")
            
    # Store results per sample
    result.df[i,"sample"] <- soi
    result.df[i,"num.loci"] <- nrow(subset_data.df)
    result.df[i,"num.match"] <- score
    result.df[i,"num.missing"] <- num_missing
    result.df[i,"prop.match"] <- prop_correct
    result.df[i, "pearson.cor"] <- calc.cor
    
    
    # then repeat for all samples
    
  }
  
  # Summarize (printout)
  print(paste0("Mean ppn of typed loci concordant per sample: ", round(mean(result.df$prop.match, na.rm = T), digits = 3)))
  print(paste0("Mean number of typed loci concordant per sample: ", round(mean(result.df$num.match, na.rm = T), digits = 3)))
  print(paste0("Mean Pearson correlation per sample: ", round(mean(result.df$pearson.cor, na.rm = T), digits = 3)))
  
  # Write out
  write.table(x = result.df, file = paste0(output_folder, "/concord_eval_", chr[c], "_comparison.txt")
              , sep = "\t", row.names = F
              )
  
  pdf(file = paste0(output_folder, "/concord_eval_", chr[c], "_hist.pdf"), width = 6.5, height = 3.5)
  hist(x = result.df$prop.match, main = "", breaks = 10
       , xlab = paste0("Per sample proportion of typed loci concordant (", chr[c], ")")
       , las = 1
       )
  #text(x = 0.67, y = 30, paste0(nrow(subset_data.df), " loci"))
  dev.off()
              
}

# Write out all info
fwrite(x = all_data.df, file = paste0(output_folder, "/all_loci_data_comparison.txt"), sep = "\t", row.names = F)

# end #
