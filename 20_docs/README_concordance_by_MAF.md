### Evaluate concordance by MAF ###
#### Setup ####
Make a working directory for this comparison:    
`mkdir eval_maf`      
Change into working dir:    
`cd eval_maf`     

Copy in the empirical data (the inds to compare to imputed), which will be used to define bins and compare to imputed.     
Copy in the imputed data, which will be used to compare.     

#### Determine bins ####
Add tags to the empirical datafile, which will be used to determine loci binned by MAF.    
```
bcftools +fill-tags <your_empirical_bcf_file>.bcf -Ob -o <your_empirical_bcf_file>_w_tags.bcf
```

Filter BCF file for MAF:     
```
# MAF bin greater than 0, less than 0.05
bcftools view -i 'MAF>0 & MAF<0.05' eval_maf/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_offspring_only_rename_w_tags.bcf -Ob -o eval_maf/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_offspring_only_rename_w_tags_maf_0.001-0.05.bcf

# MAF bin greater than 0.15 and less than 0.20
bcftools view -i 'MAF>0.15 & MAF<0.2' eval_maf/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_offspring_only_rename_w_tags.bcf -Ob -o eval_maf/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_offspring_only_rename_w_tags_maf_0.15-0.20.bcf

```

#### Comparison workflow ####
Make directories for each MAF bin output:    
e.g., `mkdir maf_to_0.05` for MAF > 0 and MAF <= 0.05.       

Conduct comparison (shown w/ 0.15 < MAF < 0.20):    
```
# Identify common loci between the MAF bin BCF and the imputed data
bcftools isec --collapse all all_inds_wgrs_and_panel_biallele_only_fi3_imputed.vcf.gz *_0.15-0.20.bcf -p maf_to_0.20

# Compress output of bcftools isec   
ls ./maf_to_0.20/*.vcf | xargs -n 1 bgzip

# Index 
ls ./maf_to_0.20/*.vcf.gz | xargs -n 1 bcftools index

# Update the run folder and target VCF.gz files of the following script and run
./01_scripts/run_bcftools_stats.sh
# ...note: following above, set
#     BCF1 as 0002.vcf.gz (imputed data present in both files); 
# and BCF2 as 0003.vcf.gz (empirical data present in both files)

# Update the target folder and run the following RScript
01_scripts/assess_bcftools_stats.R
```


