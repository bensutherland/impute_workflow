### Inspect concordance of shared loci in parents and offspring ###
Evaluate the concordance of the shared loci in the panel and wgrs data for both the parents and offspring. First, prepare the data:     
```
# Make a subfolder to keep things tidy
mkdir 05_compare/panel_vs_wgrs

# Copy in all datasets (note: filenames must be different)   
cp -l 02_input_data/parent_hd/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_parents_only_rename.bcf* 05_compare/panel_vs_wgrs/

cp -l 02_input_data/parent_ld/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_parents_only_rename.bcf* ./05_compare/panel_vs_wgrs/

cp -l 02_input_data/offspring_hd/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_offspring_only_rename.bcf* 05_compare/panel_vs_wgrs/

cp -l 02_input_data/offspring_ld/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_offspring_only_rename.bcf* 05_compare/panel_vs_wgrs/

# Merge the wgrs data back together
bcftools merge 05_compare/panel_vs_wgrs/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_parents_only_rename.bcf 05_compare/panel_vs_wgrs/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_offspring_only_rename.bcf -Ob -o 05_compare/panel_vs_wgrs/parents_and_offspring_wgrs_renamed.bcf  

bcftools index 05_compare/panel_vs_wgrs/parents_and_offspring_wgrs_renamed.bcf

# Merge the panel data back together
bcftools merge 05_compare/panel_vs_wgrs/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_parents_only_rename.bcf 05_compare/panel_vs_wgrs/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_offspring_only_rename.bcf -Ob -o 05_compare/panel_vs_wgrs/parents_and_offspring_panel_renamed.bcf

bcftools index 05_compare/panel_vs_wgrs/parents_and_offspring_panel_renamed.bcf

```

Prepare files that only contain common loci between the platforms:     
```
mkdir 05_compare/panel_vs_wgrs/isec/

# run isec to identify overlapping or private loci. Include flag '--collapse all' to consider overlap regardless of alleles.    
bcftools isec --collapse all 05_compare/panel_vs_wgrs/parents_and_offspring_wgrs_renamed.bcf 05_compare/panel_vs_wgrs/parents_and_offspring_panel_renamed.bcf -p 05_compare/panel_vs_wgrs/isec/

## Interpretation:    
# 0000.vcf = wgrs, private
# 0001.vcf = panel, private
# 0002.vcf = wgrs, shared
# 0003.vcf = panel, shared

# Save the shared output, this will be used for comparisons
bcftools view 05_compare/panel_vs_wgrs/isec/0002.vcf -Ob -o 05_compare/panel_vs_wgrs/all_inds_wgrs_shared.bcf  

bcftools index 05_compare/panel_vs_wgrs/all_inds_wgrs_shared.bcf

bcftools view 05_compare/panel_vs_wgrs/isec/0003.vcf -Ob -o 05_compare/panel_vs_wgrs/all_inds_panel_shared.bcf

bcftools index 05_compare/panel_vs_wgrs/all_inds_panel_shared.bcf

# Clean the space
rm -rf 05_compare/panel_vs_wgrs/isec
```

Make the actual comparisons
```
./01_scripts/run_bcftools_stats.sh
# ...set user variables, including target folders and whether per-site should be calculated   

```

Use the following script to analyze the outputs and generate figures:    
`01_scripts/assess_bcftools_stats.R`    

 
