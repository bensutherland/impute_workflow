### Inspect and remove Mendelian incompatibility loci ###
Mendelian incompatibilities in panel data may indicate several issues, one of which is null alleles occurring in the panel. In an imputation framework, null alleles may negatively impact the outcome of the imputation. Instructions to identify then remove loci with high frequencies of MI are given below.     

### 01. Prepare input data ###
This section expects that the panel data for parents and offspring is present in a BCF file, and that renaming may be necessary.   
The input file here will be assumed as: `02_input_data/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2.bcf`    
```
# Prepare a renaming file
bcftools query -l 02_input_data/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2.bcf > 02_input_data/samplelist.txt
# ...manually create a space-separated renaming file

# Rename the samples in the BCF file
bcftools reheader --samples 02_input_data/samplelist.txt -o 02_input_data/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_rename.bcf 02_input_data/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2.bcf

```

Use bcftools plugin Mendelian to scan for Mendelian inconsistencies, create 'bad loci' BCF file
```
# Prepare a trios text file as required by the plugin (i.e., mother1,father1,child1, 1 line per trio)  
# ...if already generated a pedigree file, can use the following shortcut
awk '{ print $3 "," $2 "," $1 }' 04_impute/pedigree.csv | grep -vE '^0' - > 06_screen_loci/trios.txt

# Use bcftools plugin Mendelian to annotate the number of Mendelian errors (MERR) in the BCF file and output
bcftools +mendelian 02_input_data/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_rename.bcf -T 06_screen_loci/trios.txt --mode a -Ob -o 06_screen_loci/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_rename_MERR.bcf

# Observe the distribution of MERR
bcftools query -f '%CHROM %POS %MERR\n' 06_screen_loci/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_rename_MERR.bcf | sort -nrk 3 | less

# Create a BCF file with problem loci
bcftools view -i 'INFO/MERR >= 4' 06_screen_loci/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_rename_MERR.bcf -Ob -o 06_screen_loci/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_rename_MERR_4.bcf

# Index
bcftools index 06_screen_loci/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_rename_MERR_4.bcf
```

Follow the main README workflow until just before imputation, then use the bad bad loci BCF file to drop bad loci from the BCF
```
# Prepare an output folder for bcftools isec
mkdir 04_impute/isec_remove_MERR

# run isec to identify loci private to the all loci data (dropping MERR)
bcftools isec 04_impute/all_inds_wgrs_and_panel_biallele.bcf 06_screen_loci/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_rename_MERR_4.bcf -p 04_impute/isec_remove_MERR/

## Interpretation:
# 0000.vcf = private to main file (no 'bad' loci)
# 0001.vcf = private to MERR file
# 0002.vcf = records from main file shared in both
# 0003.vcf = records from MERR file shared in both

# Save the private to main file vcf
bcftools view 04_impute/isec_remove_MERR/0000.vcf -Ob -o 04_impute/all_inds_wgrs_and_panel_biallele_no_MERR.bcf

# Clean up by deleting the isec folder
```

Then return to the main workflow, starting at the Imputation section.    
