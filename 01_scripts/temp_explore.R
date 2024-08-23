#### Exploring panel information, needs to be integrated into general scripts somewhere ####


parent_full_panel.vcf <- read.vcfR(file = "02_input_data/parent_ld/parent_panel.vcf")
data.df <- extract.gt(x = parent_full_panel.vcf, element = "GT")

data.df <- as.data.frame(data.df)
data.df[1:5,1:5]

for(i in 1:ncol(data.df)){
  
  print(colnames(data.df)[i])
  print(table(is.na(data.df[,i])))
  
  
}



data.df <- extract.gt(x = parent_full_panel.vcf, element = "DP")

data.df <- as.data.frame(data.df)

data.df[1:10,]

data.df <- sapply(data.df, as.numeric)
str(data.df)

colMeans(data.df)

