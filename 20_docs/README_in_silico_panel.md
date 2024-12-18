##### Optional: Simulate a panel in silico  #####
This assumes that your target file contains only the individuals targeted to be subset into panel loci, and the file has been put in `10_in_silico_panel`.     

Subsample from the target BCF file:     
```
# Obtain header
bcftools view --header-only 10_in_silico_panel/*.bcf > 10_in_silico_panel/subsample.vcf

# Subsample
bcftools view --no-header 10_in_silico_panel/*.bcf | awk '{ printf("%f\t%s\n", rand(), $0); }' | sort -t $'\t' -k1,1g | cut -f2- | head -n 1000 >> 10_in_silico_panel/subsample.vcf
# Note: solution from Pierre Lindenbaum: http://lindenb.github.io/jvarkit/DownSampleVcf.html
# Note2: here 1000 is the number of SNPs retained in the random subsample.     

# Sort and save as BCF file
bcftools sort 10_in_silico_panel/subsample.vcf -Ob -o 10_in_silico_panel/subsample_sort.bcf
bcftools index 10_in_silico_panel/subsample_sort.bcf

# This is the new in silico 'amplicon panel' file.     
```

Combine the simulated panel with parent (and grandparent) data






### OLD ###
```
# Copy in the parent data and index
bcftools index 10_impute_input/parent_wgrs/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_mi
nDP4_maxDP100_miss0.1.bcf

# Merge the subsampled offspring data with the parent data
# First see overlap
mkdir 11_impute_combine/isec_comp_sim_offsp_and_parent_wgrs

# isec to compare
bcftools isec 10_impute_input/parent_wgrs/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_min
DP4_maxDP100_miss0.1.bcf 10_impute_input/subsample_offspring_10x.bcf -p 11_impute_combine/isec_comp_sim_offsp_and_parent_wgrs/

# Interpretation
# 0000.vcf = private to wgrs
# 0001.vcf = private to simulated panel
# 0002.vcf = in wgrs in both
# 0003.vcf = in simulated panel in both

# Retain 0003.vcf
bcftools view 11_impute_combine/isec_comp_sim_offsp_and_parent_wgrs/0003.vcf -Ob -o 11_impute_combine/subsample_offspring_10x_common.bcf

# Index
bcftools index 11_impute_combine/subsample_offspring_10x_common.bcf

# Also collect a list of the loci, so this can be dropped from the 10X data later
bcftools view 11_impute_combine/subsample_offspring_10x_common.bcf | grep -vE '^#' - | awk '{ print $1 "\t" $2 }' - > 11_impute_combine/subsample_offspring_10x_common_loci_list.txt

# Combine
bcftools merge 10_impute_input/parent_wgrs/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1.bcf 11_impute_combine/subsample_offspring_10x_common.bcf -Ob -o 11_impute_combine/parent_wgrs_and_offspr_subsample.bcf

# Remove multiallelic
bcftools view --max-alleles 2 ./11_impute_combine/parent_wgrs_and_offspr_subsample.bcf -Ob -o 11_impute_combine/parent_wgrs_and_offspr_subsample_no_multiallelic.bcf
```

Next, move to the section 06. Imputation

# Note: will also want to drop the selected subsampled loci from the 10X offspring BCF that will be compared for concordance purposes, as per:
bcftools view -T ^11_impute_combine/subsample_offspring_10x_common_loci_list.txt 13_impute_compare/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1.bcf -Ob -o 13_impute_compare/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_remove_subsampled_selected_loci.bcf

# Then convert this to ai2 using the 01_scripts/bcf_to_ai2.sh script
```
