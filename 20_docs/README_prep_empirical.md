### Prepare empirical datafile to compare imputed data against ###
An empirical file will be needed to evaluate the effect of the imputation. This can be done with parents and offspring or just offspring. The offspring-only is the most appropriate comparison, given that the parents have largely full data, however, it may be of interest to see how the imputation software impacts existing complete data.     

#### 01. Prepare offspring-only file for evaluation ####
By following the main README, the file should be available here:    
`02_input_data/offspring_hd/*_offspring_only_rename.bcf`     

Note: the bcftools stats runall script `01_scripts/run_bcftools_stats.sh` will limit the comparison to only those individuals in both, as derived from the individuals that are in the empirical data (the limiting factor).     

 
