# Plotting script for manuscript only
# B. Sutherland, VIU (2024-10-04)

#### Front Matter ####
# Clean space
# rm(list=ls())

## Install and load packages
#install.packages("rstudioapi")
library(rstudioapi)
library(tidyr)
library(ggpubr)

# Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "/01_scripts", replacement = "", x = current.path)
setwd(current.path)
rm(current.path)

# Set user variable
output_folder <- "05_compare_all_loci/"

# Input materials
ai2_vs_empirical_data.FN <- "05_compare_all_loci/ai2_vs_empirical/assess_bcftools_stats_output_for_plots.RData"
fi3_vs_empirical_data.FN <- "05_compare_all_loci/fi3_vs_empirical/assess_bcftools_stats_output_for_plots.RData"

# Load AI2 plots
load(file = ai2_vs_empirical_data.FN)
# Preserve AI2 plots
hist_plot_rval_ai2         <- hist_plot_rval
hist_plot_prop_concord_ai2 <- hist_plot_prop_concord

# Load F13 plots
load(file = fi3_vs_empirical_data.FN)
# Preserve FI3 plots
hist_plot_rval_fi3         <- hist_plot_rval
hist_plot_prop_concord_fi3 <- hist_plot_prop_concord


# Create final figure
final.figure <- ggarrange(hist_plot_prop_concord_ai2  , hist_plot_rval_ai2
                          , hist_plot_prop_concord_fi3, hist_plot_rval_fi3
                          , labels = c("A", "B", "C", "D")
                          , ncol = 2, nrow = 2

)
final.figure


pdf(file = paste0(output_folder, "ai2_and_fi3_multipanel_concord.pdf"), width = 8, height = 4.5)
print(final.figure)
dev.off()



