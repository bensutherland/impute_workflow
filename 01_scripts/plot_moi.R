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
interp.FN <- "00_archive/G0923-21-VIUN_SampleInventory_V2_recd_2024-08-16.txt"

# Options
plot_targets <- "high_sig" # chr8_locus or high_sig
datatype <- "empirical" # empirical or imputed

if(datatype=="empirical"){
  
  VCF.FN    <- "plot_mois/mpileup_empirical_input_offspring_NC_047568.1.vcf.gz" # empirical
  
}else if(datatype=="imputed"){
  
  VCF.FN    <- "plot_mois/all_inds_wgrs_and_panel_biallele_only_fi3_imputed_NC_047568.1.vcf.gz" # imputed
  
}

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

if(plot_targets == "high_sig"){
  
  # High significance outliers
  moi <- c("NC_047568.1__42154072", "NC_047568.1__42154146", "NC_047568.1__42154212", "NC_047568.1__42154223" # top loci
           , "NC_047568.1__12786623", "NC_047568.1__15840611", "NC_047568.1__39651092", "NC_047568.1__35480830" # top panel loci (1)
           , "NC_047568.1__5525383",  "NC_047568.1__28802499" # top panel loci (2)
           , "NC_047568.1__19955209" # extreme outlier when no family covariate
  )
  
  
}else if(plot_targets == "chr8_locus"){ 
  
  # Loci near chr8 locus
  moi <- c("NC_047568.1__9718924", "NC_047568.1__9718930", "NC_047568.1__9719319"
           , "NC_047568.1__9719440", "NC_047568.1__9719458", "NC_047568.1__9719799"
           , "NC_047568.1__9720048")
  
  }




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
genos.df <- extract.gt(x = my_vcf, element = "GT")
dim(genos.df)
genos.df[1:5,1:5]

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
mois.df[1:5,1:5]
str(mois.df)

# Convert to df
mois.df <- as.data.frame(mois.df)

# Add indiv name
mois.df$indiv <- rownames(mois.df)
mois.df[1:5,1:5]

# Combine with interp
head(interp.df)

plot_data.df <- merge(x = interp.df, y = mois.df, by = "indiv")
dim(plot_data.df)
plot_data.df[1:5,1:10]



##### PLOT #####
# These are the available markers of interest
moi

# In the case of empirical, some of these may not actually be present in the data
present_moi <- colnames(plot_data.df)[grep(pattern = "NC_047", x = colnames(plot_data.df))]
moi <- present_moi
# # To show that its working 
# p <- ggplot(plot_data.df, aes(x = fam.x.surv, y = NC_047568.1__19955209)) + 
#   geom_boxplot(fill = "grey") + 
#   geom_jitter(height = 0.1, width = 0.25, alpha = 0.5, size = 3) + 
#   theme_classic() +
#   ylim(-0.5,2.5)
# p


# plot.list <- list()
# marker <- NULL; p <- NULL
# 
# for(i in 1:length(moi)){
#   
#   print(paste0("Plotting genotypes for marker ", moi[i]))
#   
#   marker <- moi[i]
#   marker <- sym(marker)
#   
#   p <- ggplot(plot_data.df, aes(x = fam.x.surv, y = !!marker)) + 
#     geom_boxplot(fill = "grey") + 
#     geom_jitter(height = 0.1, width = 0.25, alpha = 0.5, size = 3) + 
#     theme_classic() +
#     ylim(-0.5,2.5)
#   
#   print(p)
#   
#   plot.list[[moi[i]]] <- p
#   
# }
# 
# names(plot.list)
# print(plot.list[1])


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
    ylim(-0.5,2.5) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0)) +
    ggtitle(NULL)
  
  print(p)
  
  violin_plot.list[[moi[i]]] <- p
  
}

names(violin_plot.list)
print(violin_plot.list[1])

plot_num <- seq(1:length(violin_plot.list))

if(plot_targets=="high_sig"){
  
  final_figure <- ggarrange(violin_plot.list[[1]], violin_plot.list[[2]], violin_plot.list[[3]]
                            , violin_plot.list[[4]], violin_plot.list[[5]], violin_plot.list[[6]]
                            , violin_plot.list[[7]], violin_plot.list[[8]]
                            #, violin_plot.list[[9]]
                            #, violin_plot.list[[10]], violin_plot.list[[11]]
                            , ncol = 4, nrow = 3, labels = c("A", "B", "C", "D", "E", "F", "G", "H"
                                                             #, "I", "J", "K", "L", "M", "N"
                                                             ))
  
  pdf(file = paste0("plot_mois/multiple_loci_pheno_by_geno_violin_", plot_targets, "_", datatype, ".pdf"), width = 12, height = 9)
  print(final_figure)
  dev.off()
  
}else if(plot_targets=="chr8_locus"){
  
  final_figure <- ggarrange(violin_plot.list[[1]], violin_plot.list[[2]], violin_plot.list[[3]]
                              , violin_plot.list[[4]], violin_plot.list[[5]], violin_plot.list[[6]]
                              , violin_plot.list[[7]]
                              , ncol = 3, nrow = 3, labels = c("A", "B", "C", "D", "E", "F", "G"))
    
  pdf(file = paste0("plot_mois/multiple_loci_pheno_by_geno_violin_", plot_targets, "_", datatype, ".pdf"), width = 7, height = 7)
  print(final_figure)
  dev.off()
    
}

