### Inspect concordance of shared loci in parents and offspring ###
Evaluate concordance of shared loci in panel and wgrs data for both the parents and offspring.    

Prepare the data:     
```
# Make a subfolder to keep things tidy
mkdir 05_compare/panel_vs_wgrs

# Merge renamed parent and offspring wgrs data
bcftools merge 02_input_data/parent_hd/*_rename.bcf 02_input_data/offspring_hd/*_rename.bcf -Ob -o 05_compare/panel_vs_wgrs/parents_and_offspring_wgrs_renamed.bcf

bcftools index 05_compare/panel_vs_wgrs/parents_and_offspring_wgrs_renamed.bcf

# Merge renamed parent and offspring panel data
bcftools merge 02_input_data/parent_ld/*_rename.bcf 02_input_data/offspring_ld/*_rename.bcf -Ob -o 05_compare/panel_vs_wgrs/parents_and_offspring_panel_renamed.bcf

bcftools index 05_compare/panel_vs_wgrs/parents_and_offspring_panel_renamed.bcf

```

Prepare files with only common loci between the platforms:     
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

BCF1="all_inds_panel_shared.bcf"
BCF2="all_inds_wgrs_shared.bcf"
# note: BCF2 should be the file with the fewest individuals


```

Use the following script to analyze the outputs and generate figures:    
`01_scripts/assess_bcftools_stats.R`    

 
