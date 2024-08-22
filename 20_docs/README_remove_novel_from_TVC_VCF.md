Optional:

Remove novel variants from both panel files
```
# Remove novel variants from offspring file
# Identify hotspots (i.e., field 3 of VCF file has a marker name, and is not '.')
bcftools view 10_impute_input/offspring_panel_roslin_rehead.bcf | grep -vE '^#' - | awk '$3 != "." { print $1 "\t" $2 }' - > 10_impute_input/offspring_include_snps.txt

# Use bcftools to only keep these loci from VCF file
bcftools view --targets-file 10_impute_input/offspring_include_snps.txt ./10_impute_input/offspring_panel_roslin_rehead.bcf -Ob -o 10_impute_input/offspring_panel_roslin_rehead_hotspot_only.bcf

# Index
bcftools index 10_impute_input/offspring_panel_roslin_rehead_hotspot_only.bcf

## As above, but remove from parents file ##
bcftools view 10_impute_input/parent_panel_roslin_rehead.bcf | grep -vE '^#' - | awk '$3 != "." { print $1 "\t" $2 }' - > 10_impute_input/parent_include_snps.txt

bcftools view --targets-file 10_impute_input/parent_include_snps.txt ./10_impute_input/parent_panel_roslin_rehead.bcf -Ob -o 10_impute_input/parent_panel_roslin_rehead_hotspot_only.bcf

bcftools index 10_impute_input/parent_panel_roslin_rehead_hotspot_only.bcf
```


