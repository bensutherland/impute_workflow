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
output_folder <-"03_combine"


# Load datafiles
wgrs.vcf <- read.vcfR(file = input_wgrs_vcf.FN)
panel.vcf <- read.vcfR(file = input_panel_vcf.FN)

# Extract genotypes, then rename colnames to short form ids
wgrs_gt.df <- extract.gt(x = wgrs.vcf, element = "GT")
wgrs_gt.df[1:5,1:2]
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

# Extract genotypes, then rename colnames to short form ids
panel_gt.df <- extract.gt(x = panel.vcf, element = "GT")
panel_gt.df[1:5,1:2]
colnames(panel_gt.df)
colnames(panel_gt.df) <- gsub(pattern = "_", replacement = "-", x = colnames(panel_gt.df))
colnames(panel_gt.df) <- paste0(colnames(panel_gt.df), "_panel")
colnames(panel_gt.df)

### Join ####
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
write.table(x = all_data.df, file = paste0(output_folder, "/all_loci_data_comparison.txt"), sep = "\t", row.names = F)


