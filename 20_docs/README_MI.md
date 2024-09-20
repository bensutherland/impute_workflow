### Inspect and remove Mendelian incompatibility loci ###
Mendelian incompatibilities detected in panel data may indicate several issues, one of which is null alleles occurring in the panel. In an imputation framework, null alleles may negatively impact the outcome of the imputation. Instructions to identify then remove loci with high frequencies of MI are given below.     

### 01. Detect Mendelian incompatibilities in the panel data ###
This section expects that the panel data for parents and offspring is present in two BCF files and that they have already been renamed from IonCode to sample names.   

Prepare data:     
```
# Obtain parent and offspring LD data (note: filenames are dataset-dependent)
cp -l 02_input_data/offspring_ld/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_offspring_only_rename.bcf* 06_screen_loci/

# Merge all panel data
bcftools merge 06_screen_loci/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_parents_only_rename.bcf 06_screen_loci/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_offspring_only_rename.bcf -Ob -o 06_screen_loci/parent_and_offspring_panel.bcf

```

Detect Mendelian inconsistencies and create BCF file with 'bad loci'     
```
# Prepare a trios text file as required by the plugin (i.e., mother1,father1,child1, 1 line per trio)  
# ...if already generated a pedigree file, can use the following shortcut
awk '{ print $3 "," $2 "," $1 }' 04_impute/pedigree.csv | grep -vE '^0' - > 06_screen_loci/trios.txt

# Use bcftools plugin Mendelian to annotate the number of Mendelian errors (MERR) in the BCF file and output
bcftools +mendelian 06_screen_loci/parent_and_offspring_panel.bcf -T 06_screen_loci/trios.txt --mode a -Ob -o 06_screen_loci/parent_and_offspring_panel_MERR.bcf

# Observe the distribution of MERR
bcftools query -f '%CHROM %POS %MERR\n' 06_screen_loci/*_MERR.bcf | sort -nrk 3 | less

# Save a text file of the distribution 
bcftools query -f '%CHROM %POS %MERR\n' 06_screen_loci/*_MERR.bcf > 06_screen_loci/LD_MI_freq.txt
# ...use Rscript to plot a histogram: 01_scripts/plot_MI.R   

# Create a BCF file with problem loci
bcftools view -i 'INFO/MERR >= 4' 06_screen_loci/*_MERR.bcf -Ob -o 06_screen_loci/MERR_loci.bcf

# Index
bcftools index 06_screen_loci/MERR_loci.bcf 
```

### 02. Remove problematic loci from the analysis ###
Follow the main README workflow up to and including merging HD+LD parents with LD offspring, then use the problem loci BCF file to drop bad loci from the BCF:    
```
# Prepare an output folder for bcftools isec
mkdir 06_screen_loci/isec_rem_MERR

# Run isec to identify loci private to the all loci data (dropping MERR)
bcftools isec 04_impute/all_inds_wgrs_and_panel_biallele.bcf 06_screen_loci/MERR_loci.bcf -p 06_screen_loci/isec_rem_MERR/

## Interpretation:
# 0000.vcf = private to main file (no 'bad' loci)
# 0001.vcf = private to MERR file
# 0002.vcf = records from main file shared in both
# 0003.vcf = records from MERR file shared in both

# Save the private to main file vcf
bcftools view 06_screen_loci/isec_rem_MERR/0000.vcf -Ob -o 04_impute/all_inds_wgrs_and_panel_biallele_no_MERR.bcf

# Clean up by deleting the isec folder
rm -rf 06_screen_loci/isec_rem_MERR

# Proceed with imputation using the new *_no_MERR.bcf file
```

### 03. Optional: Detect Mendelian incompatibility loci in the wgrs data ###
Assumes that you have already separated and renamed the hd data, and it exists in `02_input_data/*_hd/`.    

```
# Merge the prepared HD data
bcftools merge 02_input_data/parent_hd/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_parents_only_rename.bcf 02_input_data/offspring_hd/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_offspring_only_rename.bcf -Ob -o 06_screen_loci/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_renamed_remerged.bcf

# See above for preparing the trios file

# Limit the trios file to only those samples that are present
bcftools query -l 06_screen_loci/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_renamed_remerged.bcf > 06_screen_loci/present_HD_samples.txt
# ...then use Rscript 01_scripts/limit_trios_file.R

# Use bcftools plugin Mendelian to annotate the number of Mendelian errors (MERR) in the BCF file and output
bcftools +mendelian 06_screen_loci/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_renamed_remerged.bcf -T 06_screen_loci/trios_limited.txt --mode a -Ob -o 06_screen_loci/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_renamed_remerged_MERR.bcf

# Observe the distribution of MERR
bcftools query -f '%CHROM %POS %MERR\n' 06_screen_loci/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_renamed_remerged_MERR.bcf > 06_screen_loci/HD_MI_freq.txt

# Use code 01_scripts/plot_MI.R to characterize

```

Then return to the main workflow, starting at the Imputation section.    
