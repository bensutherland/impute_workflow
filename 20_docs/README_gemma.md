## GWAS analysis with imputed data ##
The output of the imputation software can be used as an input to a GEMMA analysis. This is supported for FImpute3 (FI3) and AlphaImpute2 (AI2), although both are expected to be in AlphaImpute2 format (see main pipeline).       
FI3 input: `04_impute/fimpute/fi3_loci_by_inds_all_imputed_chr.txt`     
AI2 input: `05_compare/all_chr_combined.txt`      

Will also need a phenotype file (#TODO, add link) in `00_archive`.    

Use the following script to prepare the files for GEMMA, including the phenotype file, the covariate file, and the genotype file:      
`01_scripts/imputed_ai2_to_gemma.R`       

Second, run the following commands in terminal:     
```
# Change into the directory 
cd 07_GWAS/

# Construct kinship matrix
gemma -g gwas_geno.txt -p gwas_pheno.txt -gk -maf 0.05 -o gwas_all_fam

# Run GEMMA on all input files, using kinship matrix
gemma -g gwas_geno.txt -p gwas_pheno.txt -k output/gwas_all_fam.cXX.txt -n 1 -c gwas_covar.txt -maf 0.05 -lmm 4 -o gwas_all_fam_covariate

```

Third, use the following script to produce a Manhattan plot:     
`imputed_plot_gemma_results.R`     

Finally, make a directory within `07_GWAS` to hold the input and output files from your analysis, and move the files into that subfolder.   
 
