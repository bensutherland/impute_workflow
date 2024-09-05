## GWAS analysis with imputed data ##
Use the output of the imputation software, converted to ai2 format as input for this process.    
Will also need a phenotype file.    

First, use the following script to prepare the files for GEMMA, including the phenotype file, the covariate file, and the genotype file:      
`01_scripts/imputed_ai2_to_gemma.R`       

Second, run the following commands in terminal:     
```
# Change into the directory 
cd 07_impute/

# Construct kinship matrix
gemma -g gwas_geno.txt -p gwas_pheno.txt -gk -maf 0.05 -o gwas_all_fam

# Run GEMMA on all input files, using kinship matrix
gemma -g gwas_geno.txt -p gwas_pheno.txt -k output/gwas_all_fam.cXX.txt -n 1 -c gwas_covar.txt -maf 0.05 -lmm 4 -o gwas_all_fam_covariate

```

Third, use the following script to produce a Manhattan plot:     
`imputed_plot_gemma_results.R`     
- (#TODO: automate output of pval histogram)          

