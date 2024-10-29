# Compare allele frequency per locus between empirical data and imputed data (offspring only)
# B. Sutherland (2024-10-29)

# Clear space
#rm(list=ls())

# Load packages
#install.packages("rstudioapi")
#install.packages("dplyr")
#install.packages("tidyr")
#install.packages("ggplot2")
#install.packages("ggpubr")
library("tidyr")
library("rstudioapi")
library("dplyr")
library("ggplot2")
library("ggpubr")

# Set working directory to the ms_scallop_popgen repo
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path)
setwd(current.path)
rm(current.path)

## Info
# sessionInfo()
options(scipen = 99999999)

# User variables
num_to_plot <- 200


#### 01. Set up ####
# Set filenames
empirical_selected_offspring_AF.FN <- "allele_freq_comp/empirical_offspring_AF.txt"
imputed_ai2_selected_offspring_AF.FN <- "allele_freq_comp/imputed_ai2_offspring_AF.txt"
imputed_fi3_selected_offspring_AF.FN <- "allele_freq_comp/imputed_fi3_offspring_AF.txt"

## Read in data
# Empirical
emp.df <- read.delim2(file = empirical_selected_offspring_AF.FN, header = F, sep = " ")
colnames(emp.df) <- c("chr", "pos", "AF_emp")
emp.df$AF_emp <- as.numeric(emp.df$AF_emp)
emp.df$mname <- paste0(emp.df$chr, "__", emp.df$pos)
head(emp.df)

# Imputed, AI2
imp_ai2.df <- read.delim2(file = imputed_ai2_selected_offspring_AF.FN, header = F, sep = " ")
colnames(imp_ai2.df) <- c("chr", "pos", "AF_ai2")
imp_ai2.df$AF_ai2 <- as.numeric(imp_ai2.df$AF_ai2)
imp_ai2.df$mname <- paste0(imp_ai2.df$chr, "__", imp_ai2.df$pos)
head(imp_ai2.df)

# Imputed, FI3
imp_fi3.df <- read.delim2(file = imputed_fi3_selected_offspring_AF.FN, header = F, sep = " ")
colnames(imp_fi3.df) <- c("chr", "pos", "AF_fi3")
imp_fi3.df$AF_fi3 <- as.numeric(imp_fi3.df$AF_fi3)
imp_fi3.df$mname <- paste0(imp_fi3.df$chr, "__", imp_fi3.df$pos)
head(imp_fi3.df)
str(imp_fi3.df)


#### 02. Compare ####
## Empirical vs. ai2 imputed, combine
emp_v_ai2.df <- merge(x = emp.df, y = imp_ai2.df, by = "mname")
str(emp_v_ai2.df)
head(emp_v_ai2.df)
nrow(emp_v_ai2.df)

# Calculate empirical AF minus ai2 AF
emp_v_ai2.df$emp_v_ai2_AF <- emp_v_ai2.df$AF_emp - emp_v_ai2.df$AF_ai2
head(emp_v_ai2.df)

# Summary stats
summary(emp_v_ai2.df$AF_emp)
sd(emp_v_ai2.df$AF_emp)
summary(emp_v_ai2.df$AF_ai2)
sd(emp_v_ai2.df$AF_ai2)
summary(emp_v_ai2.df$emp_v_ai2_AF)
sd(emp_v_ai2.df$emp_v_ai2_AF)


# Empirical vs. fi3 imputed, combine
emp_v_fi3.df <- merge(x = emp.df, y = imp_fi3.df, by = "mname")
head(emp_v_fi3.df)
nrow(emp_v_fi3.df)

# Calculate empirical AF minus ai2 AF
emp_v_fi3.df$emp_v_fi3_AF <- emp_v_fi3.df$AF_emp - emp_v_fi3.df$AF_fi3
head(emp_v_ai2.df)

# Summary stats
summary(emp_v_fi3.df$AF_emp)
sd(emp_v_fi3.df$AF_emp)
summary(emp_v_fi3.df$AF_fi3)
sd(emp_v_fi3.df$AF_fi3)
summary(emp_v_fi3.df$emp_v_fi3_AF)
sd(emp_v_fi3.df$emp_v_fi3_AF)


# AI2 imputed
emp_v_ai2_sampled.df <- sample_n(emp_v_ai2.df, num_to_plot)

pdf(file = "allele_freq_comp/empirical_vs_alphaimpute_AF.pdf", width = 14, height = 6)
plot(y = emp_v_ai2_sampled.df$AF_emp, x = seq(1:num_to_plot), ylim = c(-0.2, 1), pch = 16, col = "blue"
     , ylab = "Allele frequency", xlab = "Index"
     , main = paste0("Empirical vs. AlphaImpute2, ", num_to_plot, " random markers")
     )
points(y = emp_v_ai2_sampled.df$AF_ai2, x = seq(1:num_to_plot), col = "orange")
abline(h = 0, lty = 3)
legend("bottomleft", legend = c("Emp.", "Imp."), col = c("blue", "orange"), pch = 19
       , cex = 0.8)
dev.off()

# Fimpute3 imputed
emp_v_fi3_sampled.df <- sample_n(emp_v_fi3.df, num_to_plot)

pdf(file = "allele_freq_comp/empirical_vs_fi3_AF.pdf", width = 14, height = 6)
plot(y = emp_v_fi3_sampled.df$AF_emp, x = seq(1:num_to_plot), ylim = c(-0.2, 1), pch = 16, col = "blue"
     , ylab = "Allele frequency", xlab = "Index"
     , main = paste0("Empirical vs. FImpute3, ", num_to_plot, " random markers")
)
points(y = emp_v_fi3_sampled.df$AF_fi3, x = seq(1:num_to_plot), col = "orange")
abline(h = 0, lty = 3)
legend("bottomleft", legend = c("Emp.", "Imp."), col = c("blue", "orange"), pch = 19
       , cex = 0.8)
dev.off()





# Things that didn't work
# pdf(file = "allele_freq_comp/empirical_vs_imputed_ai2_AF.pdf", width = 9, height = 6)
# plot(x = emp_v_ai2.df$AF_emp, y = emp_v_ai2.df$AF_ai2)
# dev.off()
# str(emp_v_fi3.df)
# test.df <- emp_v_fi3.df[1:20, ]
# 
# emp_v_fi3.contour <- ggplot(data = test.df, aes(x = AF_emp, y = AF_fi3)) +
#   stat_density_2d()
# 
# emp_v_fi3.contour


