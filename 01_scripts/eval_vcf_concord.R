# Read in two VCF files that have the same samples and check concordance
# B. Sutherland (2024-08-19)

### Front Matter ####
# Clean space
# rm(list=ls())

## Install and load packages
#install.packages("rstudioapi")
#install.packages("vcfR")
#install.packages("dartR")
library("rstudioapi")
library("vcfR")
library("dartR")
library("tidyr")

# Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "/01_scripts", replacement = "", x = current.path)
setwd(current.path)
rm(current.path)

## Info
# sessionInfo()

# Set variables
input_wgrs_vcf.FN  <- "03_combine/parent_wgrs_shared_in_both.vcf"  # VCF file of wgrs shared
input_panel_vcf.FN <- "03_combine/parent_panel_shared_in_both.vcf" # VCF file of panel shared
output_folder <-"05_compare"

input1 <- gsub(pattern = "\\.vcf", replacement = "", x = basename(input_wgrs_vcf.FN))
input2 <- gsub(pattern = "\\.vcf", replacement = "", x = basename(input_panel_vcf.FN))

# Prepare run folder
date <- format(Sys.time(), "%Y-%m-%d")
run_folder <- paste0(output_folder, "/compare_vcfs_", input1, "_vs_", input2, "_", date)
print(paste0("Making new result folder:  ", run_folder))
print(run_folder)
dir.create(run_folder)
run_folder

# Load datafiles
wgrs.vcf <- read.vcfR(file = input_wgrs_vcf.FN)
panel.vcf <- read.vcfR(file = input_panel_vcf.FN)

# Extract genotypes, file 1
wgrs_gt.df <- extract.gt(x = wgrs.vcf, element = "GT")
wgrs_gt.df[1:5,1:2]

# Rename colnames to short form ids
# TODO: this would not be needed if before running ensured the compared files had the same names
colnames.df <- as.data.frame(colnames(wgrs_gt.df))
colnames.df <- separate(data = colnames.df, col = "colnames(wgrs_gt.df)", into = c("fam", "ind", "extra")
                          , sep = "-")
head(colnames.df)
colnames(wgrs_gt.df)  <- paste0(colnames.df$fam, "-", colnames.df$ind)
rm(colnames.df)
colnames(wgrs_gt.df) <- gsub(pattern = "CH8-001", replacement = "65-8F", x = colnames(wgrs_gt.df))
colnames(wgrs_gt.df) <- gsub(pattern = "CHR8-005", replacement = "58-9F", x = colnames(wgrs_gt.df))
colnames(wgrs_gt.df)

colnames(wgrs_gt.df) <- paste0(colnames(wgrs_gt.df), "_wgrs")
colnames(wgrs_gt.df)


# Extract genotypes, file 2
# TODO: use file1 and file2 rather than wgrs and panel
panel_gt.df <- extract.gt(x = panel.vcf, element = "GT")
panel_gt.df[1:5,1:2]

# Rename colnames to short form ids
# TODO: this would not be needed if before running ensured the compared files had the same names
colnames(panel_gt.df)
colnames(panel_gt.df) <- gsub(pattern = "_", replacement = "-", x = colnames(panel_gt.df))
colnames(panel_gt.df) <- paste0(colnames(panel_gt.df), "_panel")
colnames(panel_gt.df)

### Join ####
# TODO: if use merge, then the two files would technically not need the same rownames going in
all_data.df <- cbind(wgrs_gt.df, panel_gt.df)
all_data.df <- as.data.frame(all_data.df)
dim(all_data.df)
all_data.df$mname <- rownames(all_data.df)
head(all_data.df)
all_data.df$mname <- gsub(pattern = "\\.1\\_", replacement = "\\.1__", x =  all_data.df$mname)

# Sort by colname, then put mname first
all_data.df <- all_data.df[ , order(colnames(all_data.df))]
all_data.df <- all_data.df %>% 
  select("mname", everything())
all_data.df[1:5,1:5]

# Define chr
chr <- unique(gsub(pattern = "__.*", replacement = "", x = all_data.df$mname))
chr

# Define samples
samples.vec <- unique(gsub(pattern = "_wgrs|_panel", replacement = "", x = colnames(all_data.df)))
samples.vec <- samples.vec[grep(pattern = "mname", x = samples.vec, invert = T)]
samples.vec


# Per chromosome assessment
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
  pair_compare.df <- NULL
  for(i in 1:length(samples.vec)){
    
    #print(i)
    
    # Select the sample of interest
    soi <- samples.vec[i]
    
    pair_compare.df <- subset_data.df[, c(paste0(soi, "_wgrs"), paste0(soi, "_panel"))]
    pair_compare.df <- pair_compare.df[!is.na(pair_compare.df[,1]), ]
    pair_compare.df <- pair_compare.df[!is.na(pair_compare.df[,2]), ]
    # Dropping any rows with an NA
    
    # sum up the number of identical matches between the empirical and the imputed for this sample
    score <- sum(pair_compare.df[,1] == pair_compare.df[,2])
    
    # tally the number of missing values
    num_missing <- nrow(subset_data.df) - nrow(pair_compare.df)
    
    # calculate the proportion correct for this sample by dividing the number correct by the total (with the number missing subtracted)
    prop_correct <- score / (nrow(subset_data.df) - num_missing)
    
    # # Calculate the pearson correlation between the two datatypes, only consider complete obs
    # imputed_vector   <- subset_data.df[, paste0(soi, "_imputed")]
    # empirical_vector <- subset_data.df[, paste0(soi, "_empirical")]
    # imputed_vector   <- gsub(pattern = "9", replacement = NA, x = imputed_vector)
    # empirical_vector <- gsub(pattern = "9", replacement = NA, x = empirical_vector)
    
    # imputed_vector <- as.numeric(imputed_vector)
    # empirical_vector <- as.numeric(empirical_vector)
    # 
    # calc.cor <-  cor(x = imputed_vector, y = empirical_vector
    #                  , method = "pearson", use = "pairwise.complete.obs")
    
    # Store results per sample
    result.df[i,"sample"] <- soi
    result.df[i,"num.loci"] <- nrow(subset_data.df)
    result.df[i,"num.match"] <- score
    result.df[i,"num.missing"] <- num_missing
    result.df[i,"prop.match"] <- prop_correct
    #result.df[i, "pearson.cor"] <- calc.cor
    
    
    # then repeat for all samples
    
  }
  
  # Summarize (printout)
  print(paste0("Mean ppn of typed loci concordant per sample: ", round(mean(result.df$prop.match, na.rm = T), digits = 3)))
  print(paste0("Mean number of typed loci concordant per sample: ", round(mean(result.df$num.match, na.rm = T), digits = 3)))
  #print(paste0("Mean Pearson correlation per sample: ", round(mean(result.df$pearson.cor, na.rm = T), digits = 3)))
  
  # Write out
  write.table(x = result.df, file = paste0(run_folder, "/concord_eval_", chr[c], "_comparison.txt")
              , sep = "\t", row.names = F
  )
  
  pdf(file = paste0(run_folder, "/concord_eval_", chr[c], "_hist.pdf"), width = 6.5, height = 3.5)
  hist(x = result.df$prop.match, main = "", breaks = 10
       , xlab = paste0("Per sample proportion of typed loci concordant (", chr[c], ")")
       , las = 1
  )
  #text(x = 0.67, y = 30, paste0(nrow(subset_data.df), " loci"))
  dev.off()
  
}



##### All chr together assessment #####
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
  
for(i in 1:length(samples.vec)){
    
  # Select the sample of interest
  soi <- samples.vec[i]
  
  # Select data for the sample, both VCF files
  pair_compare.df <- all_data.df[, c(paste0(soi, "_wgrs"), paste0(soi, "_panel"))]
  
  # Per locus comparison (match = TRUE)
  locus_info.df[,i+1] <- pair_compare.df[,1] == pair_compare.df[,2]
  locus_info.df[is.na(pair_compare.df[,1]), i+1] <- NA # if there is a missing record, NA the line
  locus_info.df[is.na(pair_compare.df[,2]), i+1] <- NA # if there is a missing record, NA the line
  
  # Drop any row that has missing data in either sample
  pair_compare.df <- pair_compare.df[!is.na(pair_compare.df[,1]), ]
  pair_compare.df <- pair_compare.df[!is.na(pair_compare.df[,2]), ]
  
  
  # Per individual comparison
  # sum up the number of identical matches between the empirical and the imputed for this sample
  score <- sum(pair_compare.df[,1] == pair_compare.df[,2])
  
  # tally the number of missing values
  num_missing <- nrow(all_data.df) - nrow(pair_compare.df)
  
  # calculate the proportion correct for this sample by dividing the number correct by the total (with the number missing subtracted)
  prop_correct <- score / (nrow(all_data.df) - num_missing)
  
  # Store results per sample
  result.df[i,"sample"]      <- soi
  result.df[i,"num.loci"]    <- nrow(all_data.df)
  result.df[i,"num.match"]   <- score
  result.df[i,"num.missing"] <- num_missing
  result.df[i,"prop.match"]  <- round(x = prop_correct, digits = 3)
  
}

head(result.df)
result.df


# Assess per locus
rownames(locus_info.df) <- locus_info.df$mname
locus_info.df <- locus_info.df[,2:ncol(locus_info.df)]

library("matrixStats")
locus_info.df$matches <- rowSums(x = locus_info.df, na.rm = T)
locus_info.df$NA_vals <- rowSums(is.na(locus_info.df))
locus_info.df$percent.corr <- NA

# Prop corr
for(i in 1:nrow(locus_info.df)){
  
  locus_info.df$percent.corr[i] <- locus_info.df$matches[i] / (length(samples.vec) - locus_info.df$NA_vals[i])
  
}


head(locus_info.df)



pdf(file = paste0(run_folder, "/concord_eval_per_locus_percent_corr_hist.pdf"), width = 6.5, height = 3.5)
hist(x = locus_info.df$percent.corr, main = "", breaks = 10
     , xlab = paste0("Per locus percent of complete records concordant")
     , las = 1
)
#text(x = 0.67, y = 30, paste0(nrow(subset_data.df), " loci"))
dev.off()




# Write out
write.table(x = locus_info.df, file = paste0(run_folder, "/concord_eval_per_locus_eval.txt")
            , sep = "\t", col.names = NA
)


# Write out
write.table(x = result.df, file = paste0(run_folder, "/concord_eval_all_loci_summary.txt")
            , sep = "\t", row.names = F
)
  
pdf(file = paste0(run_folder, "/concord_eval_", chr[c], "_hist.pdf"), width = 6.5, height = 3.5)
hist(x = result.df$prop.match, main = "", breaks = 10
       , xlab = paste0("Per sample proportion of typed loci concordant (", chr[c], ")")
       , las = 1
  )
  #text(x = 0.67, y = 30, paste0(nrow(subset_data.df), " loci"))
  dev.off()
  




# Write out all info
write.table(x = all_data.df, file = paste0(run_folder, "/all_loci_data_comparison.txt"), sep = "\t", row.names = F)


