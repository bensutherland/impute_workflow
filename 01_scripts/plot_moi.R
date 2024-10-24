# Extract information and genotypes for specific loci
# B. Sutherland (VIU, SBIO)
# 2024-10-24

### Front Matter ####
# Clean space
# rm(list=ls())

## Install and load packages
#install.packages("rstudioapi")
library(rstudioapi)
library(tidyr)
library(ggpubr)
library(vcfR)

# Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "/01_scripts", replacement = "", x = current.path)
setwd(current.path)
rm(current.path)

# Set user variables
VCF.FN    <- "plot_mois/all_inds_wgrs_and_panel_biallele_only_fi3_imputed_NC_047568.1.vcf.gz"
interp.FN <- "00_archive/G0923-21-VIUN_SampleInventory_V2_recd_2024-08-16.txt"


#### 01. Load data ####
# VCF file
my_vcf <- read.vcfR(file = VCF.FN)

# Interp file
interp.df <- read.table(file = interp.FN, header = T, sep = "\t")
interp.df$indiv <- gsub(pattern = "_", replacement = "-", x = interp.df$indiv)
head(interp.df)
interp.df <- separate(data = interp.df, col = "indiv", into = c("assay", "family", "rep", "id"), sep = "-", remove = F)
interp.df$fam.x.surv <- paste0(interp.df$family, "_", interp.df$survival_state)
head(interp.df)


# Drop NAs
interp.df <- interp.df[!is.na(interp.df$survival_state), ]

# Identify which outliers to plot
moi <- c("NC_047568.1__42154072", "NC_047568.1__42154146", "NC_047568.1__42154212", "NC_047568.1__42154223" # top loci
         , "NC_047568.1__12786623", "NC_047568.1__15840611", "NC_047568.1__39651092", "NC_047568.1__35480830" # top panel loci (1)
         , "NC_047568.1__5525383",  "NC_047568.1__28802499" # top panel loci (2)
)


#### 02. Extract genotypes from VCF file
# Add marker name for subsetting
head(my_vcf@fix)
my_vcf.df <- as.data.frame(my_vcf@fix)
my_vcf.df$ID <- paste0(my_vcf.df$CHROM, "__", my_vcf.df$POS)
my_vcf.df <- as.matrix(my_vcf.df)
head(my_vcf.df)
my_vcf@fix <- my_vcf.df
head(my_vcf@fix)

# Extract genotypes
genos.df <- extract.gt(x = my_vcf)
head(genos.df)

# Extract only the mois
mois.df <- genos.df[(rownames(genos.df) %in% moi), ]
dim(mois.df)
mois.df[1:5,1:5]

# Transpose
mois.df <- t(mois.df)
mois.df[1:5,1:5]

# Convert to alleles
mois.df <- gsub(pattern = "0/0", replacement = "0", x = mois.df)
mois.df <- gsub(pattern = "0/1", replacement = "1", x = mois.df)
mois.df <- gsub(pattern = "1/1", replacement = "2", x = mois.df)
mois.df[1:5,1:5]

mode(mois.df) <- "numeric"
str(mois.df)

# Convert to df
mois.df <- as.data.frame(mois.df)

# Add indiv name
mois.df$indiv <- rownames(mois.df)
mois.df[1:5,1:5]

# Combine with interp
head(interp.df)

plot_data.df <- merge(x = interp.df, y = mois.df, by = "indiv")
head(plot_data.df)

# TODO make a survival state x fam variable

##### PLOT #####
# p <- ggplot(plot_data.df, aes(x = survival_state, y = as.factor(NC_047568.1__5525383))) + 
#   geom_jitter(height = 0, width = 0.25, alpha = 0.5, size = 3) + 
#   theme_classic()
# p
# 
# p <- ggplot(plot_data.df, aes(x = survival_state, y = NC_047568.1__5525383)) + 
#   geom_boxplot(fill = "grey") +
#   theme_classic()
# p

moi

# To show that its working 
p <- ggplot(plot_data.df, aes(x = fam.x.surv, y = NC_047568.1__42154072)) + 
  geom_boxplot(fill = "grey") + 
  geom_jitter(height = 0.1, width = 0.25, alpha = 0.5, size = 3) + 
  theme_classic() +
  ylim(-0.5,2.5)
p


plot.list <- list()
marker <- NULL; p <- NULL

for(i in 1:length(moi)){
  
  print(paste0("Plotting genotypes for marker ", moi[i]))
  
  marker <- moi[i]
  marker <- sym(marker)
  
  p <- ggplot(plot_data.df, aes(x = fam.x.surv, y = !!marker)) + 
    geom_boxplot(fill = "grey") + 
    geom_jitter(height = 0.1, width = 0.25, alpha = 0.5, size = 3) + 
    theme_classic() +
    ylim(-0.5,2.5)
  
  print(p)
  
  plot.list[[moi[i]]] <- p
  
}

names(plot.list)
print(plot.list[10])


#### Violin plot option ####
violin_plot.list <- list()
marker <- NULL; p <- NULL

for(i in 1:length(moi)){
  
  print(paste0("Plotting genotypes for marker ", moi[i]))
  
  marker <- moi[i]
  marker <- sym(marker)
  
  p <- ggplot(plot_data.df, aes(x = fam.x.surv, y = !!marker)) + 
    geom_violin(fill = "grey") + 
    geom_jitter(height = 0.1, width = 0.25, alpha = 0.5, size = 3) + 
    theme_classic() +
    ylim(-0.5,2.5)
  
  print(p)
  
  violin_plot.list[[moi[i]]] <- p
  
}

names(plot.list)
print(plot.list[10])


final_figure <- ggarrange(violin_plot.list[[1]], violin_plot.list[[2]], violin_plot.list[[3]]
          , violin_plot.list[[4]], violin_plot.list[[5]], violin_plot.list[[6]]
          , violin_plot.list[[7]], violin_plot.list[[8]], violin_plot.list[[9]]
          , ncol = 3, nrow = 3, labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"))

pdf(file = "multiple_loci_pheno_by_geno_violin.pdf", width = 12, height = 12)
print(final_figure)
dev.off()








# 
# marker <- "NC_047568.1__42154072"
# marker <- sym(marker)
# 
# 
# "NC_047568.1__5525383"
# p <- ggplot(plot_data.df, aes(x = fam.x.surv, y = NC_047568.1__5525383)) + 
#   geom_boxplot(fill = "grey") + 
#   geom_jitter(height = 0.1, width = 0.25, alpha = 0.5, size = 3) + 
#   theme_classic() +  ylim(-0.5,2.5)
# p
# 
# "NC_047568.1__42154146"
# 
