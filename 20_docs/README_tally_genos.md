### Tally genotypes per sample in empirical and imputed data ###
Note: this script requires that you have rebuilt VCF files post-imputation by ai2 and fi3.    

### 01. Prepare input data
Prepare the following files to be analyzed:

```
# Combine offspring and parent data
bcftools merge 02_input_data/parent_hd/*_parents_only_rename.bcf 02_input_data/offspring_hd/*_offspring_only_rename.bcf -Oz -o 06_screen_loci/parents_and_offspring_wgrs.vcf.gz

# Obtain imputed (ai2) VCF
cp -l 04_impute/all_inds_wgrs_and_panel_biallele_only_ai2_imputed.vcf.gz ./06_screen_loci/

# Obtain imputed (fi3) VCF
cp -l 04_impute/all_inds_wgrs_and_panel_biallele_only_fi3_imputed.vcf.gz ./06_screen_loci/ 
```

### 02. Analyze proportions of different alleles ###
Use the Rscript
`01_scripts/per_sample_genotypes.R`    
... the output will be put in `06_screen_loci` and will include tables and plots.     

