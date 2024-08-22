a) Identify loci in the offspring panel file that overlap with the wgrs+panel parent file
```
# Bring offspring data to folder
cp -l 10_impute_input/offspring_panel_roslin_rehead_hotspot_only.bcf* ./11_impute_combine/

# Create isec folder to capture output
mkdir 11_impute_combine/isec_combine_parents_and_offspring/

# Use isec to compare the files (no need for --collapse here, want only common shared REF alleles)
bcftools isec 11_impute_combine/parent_wgrs_and_panel.bcf 11_impute_combine/offspring_panel_roslin_rehead_hotspot_only.bcf -p 11_impute_combine/isec_combine_parents_and_offspring/

## Interpretation:
# 0000.vcf = private to parents (wgrs+panel)
# 0001.vcf = private to offspring (panel)
# 0002.vcf = records from parents (wgrs+panel) shared in both
# 0003.vcf = records from offspring (panel) shared in both

# Save and rename 0003.vcf
cp 11_impute_combine/isec_combine_parents_and_offspring/0003.vcf 11_impute_combine/offspring_panel_roslin_rehead_hotspot_only_common_w_parents.vcf

# Compress and index
bgzip 11_impute_combine/offspring_panel_roslin_rehead_hotspot_only_common_w_parents.vcf
bcftools index 11_impute_combine/offspring_panel_roslin_rehead_hotspot_only_common_w_parents.vcf.gz
```


