### 09. Add grandparent data ###
Genotype grandparents with `wgrs_workflow`, then copy the filtered BCF to the present repository.

```
# Index the two target files before running isec
bcftools index 12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic_no_MERR.bcf
bcftools index 12_impute_impute/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP1000
0_miss0.1.bcf

# prepare an output folder for bcftools isec
mkdir 12_impute_impute/combine_all_inds_and_grandparents/

# run isec
bcftools isec 12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic_no_MERR.bcf 12_impute_impute/mpil
eup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP10000_miss0.1.bcf -p 12_impute_impute/combine_all_inds_and_grandparents/


## Interpretation:
# 0000.vcf = private to all_inds_wgrs_and_panel_no_multiallelic_no_MERR.bcf
# 0001.vcf = private to grandparents (mpileup*)
# 0002.vcf = records from all_inds shared in both
# 0003.vcf = records from grandparents shared in both

# Save out the target file to be combined
bcftools view 12_impute_impute/combine_all_inds_and_grandparents/0003.vcf -Ob -o 12_impute_impute/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP10000_miss0.1_compatible.bcf

# Index the output
bcftools index 12_impute_impute/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP10000_miss0.1_compatible.bcf

# Then can delete the isec folder to save space
```

Combine the all inds wgrs+panel with the grandparent RADseq
```
bcftools merge 12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic_no_MERR.bcf 12_impute_impute/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP10000_miss0.1_compatible.bcf -Ob -o 12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic_no_MERR_w_grandparents.bcf
```

Next, go back up to the Imputation section and run.
