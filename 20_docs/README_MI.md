#### Optional: Remove Mendelian incompatibility loci ####
Create a BCF file with parent and offspring panel-only loci:
```
# Prepare an output folder for bcftools isec
mkdir 12_impute_impute/isec_keep_only_panel_loci/

# Index
bcftools index 12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic.bcf

# run isec to identify loci shared between all loci and panel-only offspring loci
bcftools isec ./12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic.bcf 11_impute_combine/offspring
_panel_roslin_rehead_hotspot_only_common_w_parents.vcf.gz -p 12_impute_impute/isec_keep_only_panel_loci/

## Interpretation:
# 0000.vcf = private to all_inds_wgrs_and_panel_no_multiallelic.bcf
# 0001.vcf = private to offspring panel
# 0002.vcf = records from all_inds shared in both
# 0003.vcf = records from offspring panel shared in both

# Save the records from all_inds shared in both
cp -l 12_impute_impute/isec_keep_only_panel_loci/0002.vcf 12_impute_impute/all_inds_panel_only.vcf

```

Use bcftools plugin mendelian to scan for Mendelian inconsistencies
```
# Use the pedigree file that has been annotated elsewhere in this pipeline, and format as needed for the plugin
awk '{ print $3 "," $2 "," $1 }' 12_impute_impute_no_novel/pedigree_annot.csv | grep -vE '^0' - > 12_impute_impute/pedigree.csv

# Use bcftools plugin Mendelian to annotate the number of Mendelian errors (MERR) in the BCF file and output
bcftools +mendelian 12_impute_impute/all_inds_panel_only.vcf -T 12_impute_impute/pedigree.csv --mode a -Ob -o 12_impute_impute/all_inds_panel_only_annot_MERR.vcf

# Observe the distribution of MERR
bcftools query -f '%CHROM %POS %MERR\n' 12_impute_impute/all_inds_panel_only_annot_MERR.vcf | sort -nk 3 | less

# Create a BCF file with the problematic loci
bcftools view -i 'INFO/MERR >= 4' 12_impute_impute/all_inds_panel_only_annot_MERR.vcf -Ob -o 12_impute_impute/all_inds_panel_only_annot_MERR_problem_loci.bcf

# Index
bcftools index 12_impute_impute/all_inds_panel_only_annot_MERR_problem_loci.bcf
```

Remove the Mendelian inconsistencies from the wgrs+panel all individual file
```
# Prepare an output folder for bcftools isec
mkdir 12_impute_impute/isec_remove_MERR/

# run isec to identify loci private to the all loci data (dropping MERR)
bcftools isec ./12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic.bcf 11_impute_combine/offspring_panel_roslin_rehead_hotspot_only_common_w_parents.vcf.gz -p 12_impute_impute/isec_keep_only_panel_loci/

## Interpretation:
# 0000.vcf = private to all_inds_wgrs_and_panel_no_multiallelic.bcf
# 0001.vcf = private to problem loci
# 0002.vcf = records from all_inds shared in both
# 0003.vcf = records from problem loci shared in both

# Save the private records from all_inds
cp -l 12_impute_impute/isec_remove_MERR/0000.vcf 12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic_no_MERR.vcf

# Convert to BCF
bcftools view 12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic_no_MERR.vcf -Ob -o 12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic_no_MERR.bcf

# note: as clean-up, you may want to delete the isec folders, as they are large
```

