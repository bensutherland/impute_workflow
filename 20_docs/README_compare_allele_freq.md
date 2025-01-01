### Evaluate allele frequency differences between empirical and imputed genotypes ###
Generate allele frequency (AF) values for all loci common between the true (empirical) and imputed genotypes.    

NOTE: this requires that you have already conducted the standard comparison between the empirical and imputed genotypes as described in the README, and have the following files available:      
```
# FI3 imputed
05_compare/fi3_vs_empirical/all_inds_empirical_shared.bcf
05_compare/fi3_vs_empirical/all_inds_imputed_shared.bcf

# AI2 imputed
05_compare/ai2_vs_empirical/all_inds_empirical_shared.bcf
05_compare/ai2_vs_empirical/all_inds_imputed_shared.bcf
```

These files will contain all individuals in the imputed and only the present offspring in the empirical. They are already limited to contain only common loci.       


### 01. Prepare files ###
First keep only the samples from the empirical data in the imputed file:     
```
# Make a list of the samples in empirical data
bcftools query -l 05_compare/fi3_vs_empirical/all_inds_empirical_shared.bcf > 05_compare/fi3_vs_empirical/empirical_offspring_samplelist.txt

# Limit the imputed data to only contain the samples in common with the empirical
bcftools view -S 05_compare/fi3_vs_empirical/empirical_offspring_samplelist.txt 05_compare/fi3_vs_empirical/all_inds_imputed_shared.bcf -Ob -o 05_compare/fi3_vs_empirical/all_inds_imputed_shared_common_samples.bcf

```
...do the above for both fi3 and ai2 if needed.    

### 02. Obtain per-locus AF values for each file ###
Add AF values      
```
# Add AF values to all BCF files   
# FI3
bcftools +fill-tags 05_compare/fi3_vs_empirical/all_inds_empirical_shared.bcf -Ob -o 05_compare/fi3_vs_empirical/all_inds_empirical_shared_w_tags.bcf
bcftools +fill-tags 05_compare/fi3_vs_empirical/all_inds_imputed_shared_common_samples.bcf -Ob -o 05_compare/fi3_vs_empirical/all_inds_imputed_shared_common_samples_w_tags.bcf

# AI2
bcftools +fill-tags 05_compare/ai2_vs_empirical/all_inds_empirical_shared.bcf -Ob -o 05_compare/ai2_vs_empirical/all_inds_empirical_shared_w_tags.bcf
bcftools +fill-tags 05_compare/ai2_vs_empirical/all_inds_imputed_shared_common_samples.bcf -Ob -o 05_compare/ai2_vs_empirical/all_inds_imputed_shared_common_samples_w_tags.bcf

```

Extract AF values
```
# FI3
bcftools query -f '%CHROM %POS %AF\n' 05_compare/fi3_vs_empirical/all_inds_empirical_shared_w_tags.bcf > 05_compare/fi3_vs_empirical/all_inds_empirical_shared_w_tags_AF.txt
bcftools query -f '%CHROM %POS %AF\n' 05_compare/fi3_vs_empirical/all_inds_imputed_shared_common_samples_w_tags.bcf > 05_compare/fi3_vs_empirical/all_inds_imputed_shared_common_samples_w_tags_AF.txt

# AI2
bcftools query -f '%CHROM %POS %AF\n' 05_compare/ai2_vs_empirical/all_inds_empirical_shared_w_tags.bcf > 05_compare/ai2_vs_empirical/all_inds_empirical_shared_w_tags_AF.txt
bcftools query -f '%CHROM %POS %AF\n' 05_compare/ai2_vs_empirical/all_inds_imputed_shared_common_samples_w_tags.bcf > 05_compare/ai2_vs_empirical/all_inds_imputed_shared_common_samples_w_tags_AF.txt
```

### 03. Summarize, compare and plot allele frequencies per locus ###
Use the following RScript:     
`01_scripts/comp_empirical_and_imputed_AF.R`     
This will output plots of allele frequency differences for 200 random samples (set n by user), and will produce summaries of averages and differences of AFs to screen.     



