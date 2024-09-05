# Read in AlphaImpute output file
# B. Sutherland (2024-07-15)

### Front Matter ####
# Clean space
# rm(list=ls())

## Install and load packages
#install.packages("rstudioapi")
#install.packages("vcfR")
#install.packages("data.table")
#install.packages("dartR")
library("rstudioapi")
library("vcfR")
library("adegenet")
library("data.table")
library("dartR")

# Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "/01_scripts", replacement = "", x = current.path)
setwd(current.path)
rm(current.path)

## Info
# sessionInfo()

# Set variables
input.FN <- "ai2_subset.genotypes"

# Read in imputed data
imputed.mat <- fread(file = input.FN, sep = " ", header = F, stringsAsFactors = F)
dim(imputed.mat)

# Convert to df
imputed.df <- as.data.frame(imputed.mat)
imputed.df[1:5,1:5]

# Transpose
imputed_t.df <- t(imputed.df)
dim(imputed_t.df)
imputed_t.df[1:5,1:5]

# Drop ind names from matrix and use as colnames
colnames(imputed_t.df) <- imputed_t.df[1,] # use the first row as colnames
imputed_t.df <- imputed_t.df[2:nrow(imputed_t.df),]        # then drop the first row, leaving only genotypes
imputed_t.df[1:5,1:5]
imputed_t.mat <- as.matrix(imputed_t.df)
imputed_t.mat[1:5,1:5]

# Convert back to standard format
imputed_t.mat <- gsub(pattern = "\\b0\\b", replacement = "0/0", x = imputed_t.mat)
imputed_t.mat <- gsub(pattern = "\\b1\\b", replacement = "0/1", x = imputed_t.mat)
imputed_t.mat <- gsub(pattern = "\\b2\\b", replacement = "1/1", x = imputed_t.mat)

# Read in the contributing VCF file
input_vcf.FN <- "12_impute_impute/all_inds_wgrs_and_panel_NC_047567_1.vcf"


## Import VCF for mnames
input.vcf <- vcfR::read.vcfR(file = input_vcf.FN)
input.vcf

# # Obtain genotypes too
# test.df <- extract.gt(x = input.vcf, element = "GT")
# test.df[1:5,1:5]
# test.df[1:5, 20:25]
# colnames(test.df)

# Get mnames
dim(input.vcf@gt)
input.vcf@gt[1:5,1:5]
dim(input.vcf@gt)
input.vcf@fix[1:5,]

geno_info.df <- input.vcf@fix
geno_info.df <- as.data.frame(geno_info.df)

# Generate an mname based on the VCF file
mname <- paste0(geno_info.df$CHROM, "__", geno_info.df$POS, "__", geno_info.df$ID)
length(mname)

imputed_t.mat[1:5,1:5]
rownames(imputed_t.mat) <- mname

imputed_t.mat[1:5,1:5]
mnames.df <- rownames(imputed_t.mat)
mnames.df <- as.data.frame(mnames.df)

imputed_t.df <- as.data.frame(imputed_t.mat)

imputed_full.df <- cbind(mnames.df, imputed_t.df)

fwrite(x = imputed_full.df, file = "12_impute_impute/genos_imputed_converted.txt", sep = "\t", quote = F)

# Transpose again
data.mat <- t(imputed_t.mat) 
data.mat[1:5,1:5]

# Remove any periods in the header
colnames(data.mat) <- gsub(pattern = "\\.", replacement = "_", x = colnames(data.mat))

# Convert to genind file
data.gid <- df2genind(X = data.mat, sep = "/")

# Convert genind to genlight
data.gl <- gi2gl(gi = data.gid, parallel = T)

pca1 <- glPca(x = data.gl, nf = 4, parallel = T, loadings = F)

# Set working directory to the ms_scallop_popgen repo
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "/01_scripts", replacement = "", x = current.path)
setwd(current.path)
rm(current.path)

pca_from_genind(data = data.gid, PCs_ret = 4, plot_eigen = F, plot_allele_loadings = F
                , plot_ellipse = F, retain_pca_obj = T
                )
#pca1.bck <- pca1

# pop.id <- indNames(data.gid)
# pop.id[grep(pattern = "ASY2", x = pop.id, invert = T)] <- "F0"
# pop.id[grep(pattern = "_114_", x = pop.id)] <- "F114"
# pop.id[grep(pattern = "_115_", x = pop.id)] <- "F115"
# pop.id[grep(pattern = "_116_", x = pop.id)] <- "F116"
# pop.id[grep(pattern = "_117_", x = pop.id)] <- "F117"
# pop.id[grep(pattern = "55-41F")]
# pop.id

pca_scores.df <- pca1$scores
pca_scores.df <- as.data.frame(pca_scores.df)

pca_scores.df$pop <- rownames(pca_scores.df)
head(pca_scores.df)
pca_scores.df$pop[grep(pattern = "ASY2", x = pca_scores.df$pop, invert = T)] <- "F0"
pca_scores.df$pop[grep(pattern = "_114_", x = pca_scores.df$pop)] <- "F114"
pca_scores.df$pop[grep(pattern = "_115_", x = pca_scores.df$pop)] <- "F115"
pca_scores.df$pop[grep(pattern = "_116_", x = pca_scores.df$pop)] <- "F116"
pca_scores.df$pop[grep(pattern = "_117_", x = pca_scores.df$pop)] <- "F117"
pca_scores.df$pop

# Plot
p <- ggplot(data = pca_scores.df, aes(x = PC1, y = PC2, colour=pop))
p <- p + geom_point(size = 2)
p

pdf(file = "12_impute_impute/subset_chr_pca.pdf", width = 6.5, height = 5.5)
print(p)
dev.off()

write.table(x = pca_scores.df, file = "12_impute_impute/subset_chr_pca.txt", sep = "\t", quote = F
            , col.names = NA, row.names = T
            )

# Next: gemma preparation 
