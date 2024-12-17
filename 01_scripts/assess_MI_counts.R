# Assess MI results
# B. Sutherland (2024-12-17)

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

# Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "/01_scripts", replacement = "", x = current.path)
setwd(current.path)
rm(current.path)

# Set variables
#MI.FN     <- "06_screen_loci/parents_and_offspring_wgrs_MI_count.txt"
MI.FN     <- "06_screen_loci/parents_and_offspring_panel_MI_count.txt"

# Read in results
MI.df <- read.table(file = MI.FN, sep = "\t", header = F)
head(MI.df)
colnames(MI.df) <- c("nOK", "nBad", "nSkipped", "trio")
head(MI.df)

# Separate into families
MI.df$fam <- NA
MI.df$fam[grep(pattern = "-114-", x = MI.df$trio)] <- "F114"
MI.df$fam[grep(pattern = "-115-", x = MI.df$trio)] <- "F115"
MI.df$fam[grep(pattern = "-116-", x = MI.df$trio)] <- "F116"
MI.df$fam[grep(pattern = "-117-", x = MI.df$trio)] <- "F117"
table(MI.df$fam)

# Add metric columns
MI.df$prop.scored.bad <- MI.df$nBad / (MI.df$nOK + MI.df$nBad)
hist(MI.df$prop.scored.bad)

MI.df$prop.skipped <- MI.df$nSkipped / (MI.df$nOK + MI.df$nBad + MI.df$nSkipped)
hist(MI.df$prop.skipped)

#plot(MI.df$prop.scored.bad ~ MI.df$prop.skipped, las = 1, col = as.factor(MI.df$fam))

# Plot
p <- ggplot(MI.df, aes(x=prop.skipped, y=prop.scored.bad, colour=fam)) + 
       geom_point() + 
        xlab("Proportion missing data") + 
        ylab("Proportion of loci MI")
p

pdf(file = gsub(pattern = ".txt", replacement = "_scatterplot.pdf", x = MI.FN), width = 5.5, height = 4)
print(p)
dev.off()

head(MI.df)

# Summarize by family, average MI proportion
MI.df %>% 
  group_by(fam) %>%
  summarize(mean.prop.bad = mean(prop.scored.bad))

# Summarize by family, median MI proportion
MI.df %>% 
  group_by(fam) %>%
  summarize(median.prop.bad = median(prop.scored.bad))

# Summarize by family, sd MI proportion
MI.df %>% 
  group_by(fam) %>%
  summarize(sd.prop.bad = sd(prop.scored.bad))

# Summarize by family, avg. prop skipped
MI.df %>% 
  group_by(fam) %>%
  summarize(mean.prop.skipped = mean(prop.skipped))

# Summarize by family, avg. ok
MI.df %>% 
  group_by(fam) %>%
  summarize(mean.nOK = mean(nOK))

# Summarize by family, avg. bad
MI.df %>% 
  group_by(fam) %>%
  summarize(mean.nBad = mean(nBad))

# Summarize by family, avg. skipped
MI.df %>% 
  group_by(fam) %>%
  summarize(mean.nskipped = mean(nSkipped))

save.image(file = gsub(pattern = ".txt", replacement = "_output.RData", x = MI.FN))
