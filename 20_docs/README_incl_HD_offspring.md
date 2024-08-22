#### Optional: Include some high density offspring to support the imputation ####
Copy in the all inds wgrs and panel BCF file, as well as the all offspring wgrs 10x data into `12_impute
_impute`.
Prepare files for combining:
```
# Remove four offspring from each family from the imputation target file
bcftools view 12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic_no_MERR.bcf --samples ^ASY2_114_R
1_1,ASY2_114_R1_2_ReAMP,ASY2_114_R1_3_ReAMP,ASY2_114_R1_4_ReAMP,ASY2_115_R1_1_ReAMP,ASY2_115_R1_2_ReAMP,ASY2_115_R1_3_ReAMP,ASY2_115_R1_4_ReAMP,ASY2_116_R1_1_ReAMP,ASY2_116_R1_2_ReAMP,ASY2_116_R2_1_ReAMP,ASY2_116_R2_2_ReAMP,ASY2_117_R1_1_ReAMP,ASY2_117_R1_2_ReAMP,ASY2_117_R1_3_ReAMP,ASY2_117_R1_4_ReAMP -Ob -o 12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic_no_MERR_rem_16_inds.bcf

# Index
bcftools index 12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic_no_MERR_rem_16_inds.bcf


# Isolate these same individuals from the 10X wgrs offspring file
bcftools view 12_impute_impute/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1.bcf --samples ASY2-114-R1-1-43820254_S191_L002_R1.fastq.gz,ASY2-114-R1-2-43820255_S192_L002_R1.fastq.gz,ASY2-114-R1-3-43820256_S193_L002_R1.fastq.gz,ASY2-114-R1-4-43820257_S194_L002_R1.fastq.gz,ASY2-115-R1-1-43820314_S246_L002_R1.fastq.gz,ASY2-115-R1-2-43820315_S247_L002_R1.fastq.gz,ASY2-115-R1-3-43820316_S248_L002_R1.fastq.gz,ASY2-115-R1-4-43820317_S249_L002_R1.fastq.gz,ASY2-116-R1-1-43820371_S31_L002_R1.fastq.gz,ASY2-116-R1-2-43820372_S32_L002_R1.fastq.gz,ASY2-116-R2-1-43820379_S34_L002_R1.fastq.gz,ASY2-116-R2-2-43820380_S35_L002_R1.fastq.gz,ASY2-117-R1-1-43820428_S63_L002_R1.fastq.gz,ASY2-117-R1-2-43820429_S64_L002_R1.fastq.gz,ASY2-117-R1-3-43820430_S65_L002_R1.fastq.gz,ASY2-117-R1-4-43820431_S66_L002_R1.fastq.gz -Ob -o 12_impute_impute/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_select_16_inds.bcf

# Index
bcftools index 12_impute_impute/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_select_16_inds.bcf

# Compare the two prepared files via isec
mkdir 12_impute_impute/isec_impute_target_vs_10X/

bcftools isec 12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic_no_MERR_rem_16_inds.bcf 12_impute_impute/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_select_16_inds.bcf -p 12_impute_impute/isec_impute_target_vs_10X/

## Interpretation:
# 0000.vcf = private to impute target file
# 0001.vcf = private to 10x file
# 0002.vcf = records from impute target file shared in both
# 0003.vcf = records from 10x file shared in both

# Use 0003.vcf, but convert to bcf first, and save out
bcftools view 12_impute_impute/isec_impute_target_vs_10X/0003.vcf -Ob -o 12_impute_impute/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_select_16_inds_w_compat_loci.bcf

bcftools index 12_impute_impute/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_select_16_inds_w_compat_loci.bcf

# Clean up space by deleting the isec folder

# Combine
bcftools merge 12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic_no_MERR_rem_16_inds.bcf 12_impute_impute/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_select_16_inds_w_compat_loci.bcf -Ob -o 12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic_no_MERR_add_16_HD_offspring.bcf
```
