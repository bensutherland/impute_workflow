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
output_folder <- "05_compare/"

# Input materials
panel_vs_wgrs_data.FN    <- "05_compare/panel_vs_wgrs/assess_bcftools_stats_output_for_plots.RData"
ai2_vs_empirical_data.FN <- "05_compare/ai2_vs_empirical/assess_bcftools_stats_output_for_plots.RData"
fi3_vs_empirical_data.FN <- "05_compare/fi3_vs_empirical/assess_bcftools_stats_output_for_plots.RData"

## Load panel vs wgrs data
load(file = panel_vs_wgrs_data.FN)

# Preserve panel_vs_wgrs plots
hist_plot_prop_panel_vs_wgrs   <- hist_plot_prop_concord
hist_plot_rval_panel_vs_wgrs   <- hist_plot_rval
psd_plot_panel_vs_wgrs         <- psd_plot

# Fix xlim
hist_plot_rval_panel_vs_wgrs   <- hist_plot_rval_panel_vs_wgrs + xlim(0.4, 1) # fix


## Load AI2 plots
load(file = ai2_vs_empirical_data.FN)

# Preserve AI2 plots
hist_plot_prop_concord_ai2  <- hist_plot_prop_concord
hist_plot_rval_ai2          <- hist_plot_rval
psd_plot_ai2 <- psd_plot


## Load F13 plots
load(file = fi3_vs_empirical_data.FN)

# Preserve FI3 plots
hist_plot_prop_concord_fi3  <- hist_plot_prop_concord
hist_plot_rval_fi3          <- hist_plot_rval
psd_plot_fi3                <- psd_plot


#### Create final figure, no per-locus vals ####
final.figure <- ggarrange(
  hist_plot_prop_panel_vs_wgrs, hist_plot_rval_panel_vs_wgrs
  , hist_plot_prop_concord_ai2  , hist_plot_rval_ai2
  , hist_plot_prop_concord_fi3  , hist_plot_rval_fi3
  , labels = c("A", "B", "C", "D", "E", "F")
  , ncol = 2, nrow = 3
  
)
final.figure


pdf(file = paste0(output_folder, "panel_vs_wgrs_ai2_and_fi3_imputed_vs_empirical_multipanel.pdf"), width = 7.5, height = 4.5)
print(final.figure)
dev.off()






# #### Create final figure (nine panel) ####
# final.figure <- ggarrange(
#     hist_plot_prop_panel_vs_wgrs, hist_plot_rval_panel_vs_wgrs, psd_plot_panel_vs_wgrs
#   , hist_plot_prop_concord_ai2  , hist_plot_rval_ai2, psd_plot_ai2
#   , hist_plot_prop_concord_fi3  , hist_plot_rval_fi3, psd_plot_fi3
#   , labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I")
#   , ncol = 3, nrow = 3
# 
# )
# final.figure
# 
# 
# pdf(file = paste0(output_folder, "panel_vs_wgrs_ai2_and_fi3_imputed_vs_empirical_multipanel.pdf"), width = 9, height = 5.5)
# print(final.figure)
# dev.off()


# # Create final figure (four panel option)
# final.figure <- ggarrange(hist_plot_prop_concord_ai2  , hist_plot_rval_ai2
#                           , hist_plot_prop_concord_fi3, hist_plot_rval_fi3
#                           , labels = c("A", "B", "C", "D")
#                           , ncol = 2, nrow = 2
# 
# )
# final.figure
# 
# 
# pdf(file = paste0(output_folder, "ai2_and_fi3_multipanel_concord.pdf"), width = 8, height = 4.5)
# print(final.figure)
# dev.off()
