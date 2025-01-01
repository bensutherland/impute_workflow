# impute_workflow
Imputation workflow for working with amplicon panel and wgrs data

### 00. Getting started ###
Clone the present repo, all commands will occur in the repo unless indicated.    

Requires the following inputs put in `02_input_data`:      
- high-density (wgrs) BCF file
- low-density (panel) BCF file

The default pipeline assumes that the high-density and low-density data were generated independently, and that within each are parents and offspring genotyped together. Some support of different inputs is provided (#TODO, link).     

Note: no duplicate individuals (i.e., tech replicates) should be present.    

#### Special Functions ####
- Screen for Mendelian incompatibilities [here](20_docs/README_MI.md)     
- Compare between platform concordance [here](20_docs/README_compare_shared_loci_HD_LD.md)    
- Screen for unexpected parental monomorphism per family [here](20_docs/README_unexpected_monomorphism.md)      
- Tally genotypes per sample in empirical and imputed data [here](20_docs/README_tally_genos.md)    
- Compare allele freq. per locus in empirical and imputed data [here](20_docs/README_compare_allele_freq.md)      


### 01. Prepare input data ### 
#### wgrs data (HD data) ####
##### Separate and rename parents and offspring #####   
Separate and rename parents:      
```
# Create a parent samplefile, assumes non-parents all have 'ASY2' string in name
bcftools query -l 02_input_data/<wgrs_data>.bcf | grep -vE '^ASY2' - > 02_input_data/parent_hd_samplelist.txt
# ...note: if there are any other parents you want to remove from the dataset, delete them from this select list here.     

# Select only the parents from the dataset
bcftools view -S 02_input_data/parent_hd_samplelist.txt 02_input_data/<wgrs_data>.bcf -Ob -o 02_input_data/parent_hd/<wgrs_data>_parents_only.bcf

# Rename 
# ...manually add desired samplenames to 02_input_data/parent_hd_samplelist.txt in space-separated format
bcftools reheader --samples 02_input_data/parent_hd_samplelist.txt -o 02_input_data/parent_hd/<wgrs_data>_parents_only_rename.bcf 02_input_data/parent_hd/<wgrs_data>_parents_only.bcf

# Index
bcftools index 02_input_data/parent_hd/<wgrs_data>_parents_only_rename.bcf

# Clean space
rm 02_input_data/parent_hd/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_parents_only.bcf
```

Separate and rename offspring:      
```
# Create an offspring samplefile, assumes non-parents all have 'ASY2' string in name
bcftools query -l 02_input_data/<wgrs_data>.bcf | grep -E '^ASY2' - > 02_input_data/offspring_hd_samplelist.txt
# ...note: if there are any other offspring you want to remove from the dataset, delete them from this select list here

# Select only the offspring from the dataset
bcftools view -S 02_input_data/offspring_hd_samplelist.txt 02_input_data/<wgrs_data>.bcf -Ob -o 02_input_data/offspring_hd/<wgrs_data>_offspring_only.bcf

# Rename
# ...here demonstrated in automated method from existing samplenames
awk -F"-" '{ print $0 " " $1"-"$2"-"$3"-"$4 }' 02_input_data/offspring_hd_samplelist.txt > 02_input_data/offspring_hd_samplelist_for_rename.txt

bcftools reheader --samples 02_input_data/offspring_hd_samplelist_for_rename.txt -o 02_input_data/offspring_hd/<wgrs_data>_offspring_only_rename.bcf 02_input_data/offspring_hd/<wgrs_data>_offspring_only.bcf

# Index
bcftools index 02_input_data/offspring_hd/<wgrs_data>_offspring_only_rename.bcf

# Clean space
rm 02_input_data/offspring_hd/<wgrs_data>_offspring_only.bcf
```

#### panel data (LD data) ####
##### Separate and rename parents and offspring #####   
Separate and rename parents:      
```
# Create a parent samplefile (in this example, identified by string 'rawlib')
bcftools query -l 02_input_data/<panel_data>.bcf | grep 'rawlib' - > 02_input_data/parent_ld_samplelist.txt
# ...note: if there are any parents you want to remove, delete from this samplelist

# Select only the parents from the dataset
bcftools view -S 02_input_data/parent_ld_samplelist.txt 02_input_data/<panel_data>.bcf -Ob -o 02_input_data/parent_ld/<panel_data>_parents_only.bcf

# Rename
# ...manually add desired samplenames to 02_input_data/parent_ld_samplelist.txt in space-separated format
bcftools reheader --samples 02_input_data/parent_ld_samplelist.txt -o 02_input_data/parent_ld/<panel_data>_parents_only_rename.bcf 02_input_data/parent_ld/<panel_data>_parents_only.bcf

# Index
bcftools index 02_input_data/parent_ld/<panel_data>_parents_only_rename.bcf
```

Separate and rename offspring:    
```
# Create an offspring samplefile (in this example, identified by absence of 'rawlib')
bcftools query -l 02_input_data/<panel_data>.bcf | grep -vE 'rawlib' - > 02_input_data/offspring_ld_samplelist.txt
# ...note: if there are any offspring you want to remove, delete from this samplelist

# Select only the offspring from the dataset
bcftools view -S 02_input_data/offspring_ld_samplelist.txt 02_input_data/<panel_data>.bcf -Ob -o 02_input_data/offspring_ld/<panel_data>_offspring_only.bcf

# Rename
# ...manually create 02_input_data/offspring_ld_samplelist_for_rename.txt
bcftools reheader --samples 02_input_data/offspring_ld_samplelist_for_rename.txt -o 02_input_data/offspring_ld/<panel_data>_offspring_only_rename.bcf 02_input_data/offspring_ld/<panel_data>_offspring_only.bcf

# Index
bcftools index 02_input_data/offspring_ld/<panel_data>_offspring_only_rename.bcf
```

Before leaving this section, create a useful file that can be used later, for example in plotting, that indicates what SNPs were from the panel originally:     
```
# Collect the names (i.e., 'chr__pos') of the panel variants
bcftools view 02_input_data/offspring_ld/<panel_data>_offspring_only_rename.bcf | grep -vE '^#' - | awk '{ print $1 "__" $2 }' - > 07_GWAS/denovo_snp_ids.txt
```


### 02. Exclude panel loci from parent wgrs data ###
Before merging parent wgrs and panel data, all overlapping loci between the two will be dropped from the wgrs dataset. Note: this will also give substrate for direct comparison of overlapping loci (see next section).     

Separate loci into private or overlapping:     
```
# prepare an output folder for bcftools isec
mkdir 03_combine/isec_wgrs_and_panel

# run isec to identify overlapping or private loci. Include flag '--collapse all' to consider overlap regardless of alleles.    
bcftools isec --collapse all 02_input_data/parent_hd/<wgrs_data>_parents_only_rename.bcf 02_input_data/parent_ld/<panel_data>_parents_only_rename.bcf -p 03_combine/isec_wgrs_and_panel/

## Interpretation:    
# 0000.vcf = wgrs, private
# 0001.vcf = panel, private
# 0002.vcf = wgrs, shared
# 0003.vcf = panel, shared

# Save the wgrs, private VCF file as BCF file
bcftools view 03_combine/isec_wgrs_and_panel/0000.vcf -Ob -o 03_combine/parent_wgrs_only.bcf

# Save the overlapped loci as VCF files for concordance evaluation (next section)
cp -l 03_combine/isec_wgrs_and_panel/0002.vcf 05_compare/parent_wgrs_shared_in_both.vcf
cp -l 03_combine/isec_wgrs_and_panel/0003.vcf 05_compare/parent_panel_shared_in_both.vcf

# Delete the isec folder to save space

```


### 03. Concatenate parent panel loci into parent wgrs only data ###
Prepare the parent wgrs-only file to be combined:       
```
# Create a sorted samplelist text file
bcftools query -l 03_combine/parent_wgrs_only.bcf | sort > 03_combine/parent_wgrs_samplenames_sorted.txt

# Sort in the BCF file
bcftools view -S 03_combine/parent_wgrs_samplenames_sorted.txt 03_combine/parent_wgrs_only.bcf -Ob -o 03_combine/parent_wgrs_only_sorted.bcf

# Index 
bcftools index 03_combine/parent_wgrs_only_sorted.bcf
```

Prepare the parent panel data to be combined:       
```
# Copy the parent panel data into the combined folder
cp -l 02_input_data/parent_ld/<panel_data>_parents_only_rename.bcf 03_combine/

# Create sorted text file of names 
bcftools query -l 03_combine/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_parents_only_rename.bcf | sort > 03_combine/parent_panel_samplenames_sorted.txt

# Sort in the BCF file
bcftools view -S 03_combine/parent_panel_samplenames_sorted.txt 03_combine/<panel_data>_parents_only_rename.bcf -Ob -o 03_combine/<panel_data>_parents_only_rename_sorted.bcf

# index
bcftools index 03_combine/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_parents_only_rename_sorted.bcf
```

Combine the parent data with bcftools concat
```
bcftools concat --allow-overlaps 03_combine/parent_wgrs_only_sorted.bcf 03_combine/<panel_data>_parents_only_rename_sorted.bcf -Ob -o 03_combine/parent_wgrs_and_panel.bcf

# Index
bcftools index 03_combine/parent_wgrs_and_panel.bcf
```


### 04. Merge parent wgrs and panel data with offspring panel data ###
Combine the wgrs+panel parent data with the panel offspring data
```
# Merge parent wgrs+panel with offspring panel
bcftools merge 03_combine/parent_wgrs_and_panel.bcf 02_input_data/offspring_ld/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_offspring_only_rename.bcf -Ob -o 03_combine/all_inds_wgrs_and_panel.bcf

# Remove multiallelic if present
bcftools view --max-alleles 2 03_combine/all_inds_wgrs_and_panel.bcf -Ob -o 03_combine/all_inds_wgrs_and_panel_biallele.bcf

# Copy to imputation folder
cp -l 03_combine/all_inds_wgrs_and_panel_biallele.bcf ./04_impute/

# Index 
bcftools index 04_impute/all_inds_wgrs_and_panel_biallele.bcf
```


### 05. Imputation ###
The data is now all in a single BCF file and is ready for the imputation process.     

Prepare pedigree file:      
```
bcftools query -l 04_impute/all_inds_wgrs_and_panel_biallele.bcf > 04_impute/pedigree.csv

# Annotate the above file as follows:    
# <indiv> <sire> <dam>     
# where if there is no sire or dam, put 0
# save as space-delimited    
```

#### Imputation with AlphaImpute2: ####   
Prepare data:    
```
# Convert from BCF file to AlphaImpute2 format
./01_scripts/bcf_to_ai2.sh
# produces: 04_impute/<input_filename>_ai2.txt
# in format:  mname \t ind1 \t ind2 \t (...)
#             NC_047559.1 2945 \t 0 \t 0 \t 1 \t 9 (...)
# where 0, 1, 2 are the number of alt alleles, and 9 is missing data  

# Separate ai2 matrix into individual chr. Set input filename, and string to identify chromosomes.
01_scripts/prep_geno_matrix_for_ai2.R   
# produces: 04_impute/ai2_input_<NC_047559.1>.txt, one file per chr 
```

Impute:    
```
# initialize the conda environment
conda activate ai2

# Run AlphaImpute2 on chromosome-separated datafiles
01_scripts/run_ai2.sh
# produces: 04_impute/ai2_input_<NC_047559.1>.genotypes and .haplotypes

# Transpose chr-sep .genotypes files, and prepare to combine back together
01_scripts/impute_rebuild_chr_lightweight.R
# produces: 05_compare/*.genotypes_transposed_to_combine.txt

# Combine imputed, transposed, chr-separated files back together, then add marker names back in, based on the input ai2 file (before chromosome separation) 
#   note: before running, update the variable for your input original ai2 file
01_scripts/combine_transposed_ai2_output_and_mnames.sh
# produces: 05_compare/all_chr_combined.txt

```

Rebuild VCF file from ai2 output:    
```
# Need to create a VCF file with only the loci present in the imputed output
# First, identify which regions were output by ai2
awk '{ print $1 "\t" $2  }' 05_compare/all_chr_combined.txt | grep -vE '^mname' - > 05_compare/ai2_imputed_regions.txt

# Second, subset the pre-imputed VCF file to only include these selected regions
bcftools view --regions-file 05_compare/ai2_imputed_regions.txt 04_impute/all_inds_wgrs_and_panel_biallele.bcf -Ov -o 04_impute/all_inds_wgrs_and_panel_biallele_only_ai2_imputed_regions.vcf

# Using ai2 output file, convert allele dosage (0,1,2) format to standard VCF format (0/0, 0/1, 1/1)
./01_scripts/ai2_to_vcf_format_step_1.sh
# ...output will be the input identifier, with "_converted.txt" as suffix

# Use the following Rscript to read in the VCF and replace the genotypes in the pre-imputed with the imputed data
./01_scripts/ai2_to_VCF.R

# Convert the output to BCF format
bcftools view 04_impute/<your_output_file>.vcf.gz -Ob -o 04_impute/<your_output_file>_ai2_imputed.bcf

# Index the above BCF file
```

Jump to [Evaluate concordance](#07-evaluate-imputation).      



#### Imputation with FImpute3: ####   
The FImpute3 workflow will start from the pre-imputed, pre-chromosome separated AI2 file generated from the BCF file above.     
 
Further prepare data:    
```
# Convert from ai2 to FImpute3 format
./01_scripts/ai2_to_fimpute3.R
# ...this creates a pedigree file, a map file per chr, and parts needed for genotypes file per chr

# *Note*: this will assume females and males sample names end with F or M, respectively, otherwise will default to male. If there are parents without this name format, edit the pedigree file that was produced in 04_impute/fimpute/pedigree.csv 
# Note: if you are editing, ensure that the order is sire then dam for fimpute3, else an error will occur

# Edit and combine components needed for per-chr genotypes files
./01_scripts/prep_fi3_indiv_chr.sh

# Prepare a control file and output run directory per chr 
./01_scripts/prep_fi3_control_files_and_output_dir.sh

```

Impute:    
```
# Run FImpute3 iteratively
./01_scripts/run_fi3_iteratively.sh

# Convert the FImpute3 output back to ai2 format to prepare for conversion to VCF
./01_scripts/fimpute3_to_ai2.R
```

Rebuild VCF file from fi3 output:    
```
# Create a VCF file containing only those loci present after imputation
# First, identify which regions were output by imputation program
awk '{ print $1 "\t" $2 }' 04_impute/fimpute/fi3_loci_by_inds_all_imputed_chr.txt | grep -vE '^mname' - > 05_compare/fi3_imputed_regions.txt

# Second, subset the pre-imputed VCF file to only include these selected regions
bcftools view --regions-file 05_compare/fi3_imputed_regions.txt 04_impute/all_inds_wgrs_and_panel_biallele.bcf -Ov -o 04_impute/all_inds_wgrs_and_panel_biallele_only_fi3_imputed_regions.vcf

# Using the output imputation file, convert from allele dosage (0,1,2) format to standard VCF format (0/0, 0/1, 1/1)
# note: set input folder and input ai2-formatted file
./01_scripts/ai2_to_vcf_format_step_1.sh
# ...output will be the input identifier, with "_converted.txt" as suffix

# Use the following Rscript to read in the VCF and replace the genotypes in the pre-imputed with the imputed data
./01_scripts/ai2_to_VCF.R
# ...output will be the input vcf (w/ regions) identifier, with "_fi3|_ai2_imputed.vcf.gz" 

# Convert the output to BCF format
bcftools view 04_impute/<your_output_file>.vcf.gz -Ob -o 04_impute/all_inds_wgrs_and_panel_biallele_only_fi3_imputed.bcf

# Index the above BCF file
```

Jump to [Evaluate concordance](#07-evaluate-imputation).      


### 07. Evaluate imputation ###
This section will describe how to compare BCF/VCF files.      

See [here](20_docs/README_prep_empirical.md) for details on how to obtain a wgrs empirical file.     
See [here](20_docs/README_compare_shared_loci_HD_LD.md) for details on how to do a between platform comparison (low density vs. high density shared loci).     

Prepare data for comparison:          
note: use a name below that fits your comparison type, for example `ai2_vs_empirical`, `fi3_vs_empirical`, `wgrs_vs_panel`      
```
# note: this example will use ai2 vs empirical, but other approaches follow the same

# Make a subfolder to keep things tidy
mkdir 05_compare/ai2_vs_empirical

# Copy in both files, for example: 
# a prepared file that has all wgrs genotypes (samples renamed)
# 04_impute/all_inds_wgrs_and_panel_biallele_only_ai2_imputed.bcf (imputed data)
# 02_input_data/offspring_hd/*_offspring_only_rename.bcf (empirical data)

# Prepare an isec folder
mkdir 05_compare/ai2_vs_empirical/isec/

# Run isec
bcftools isec --collapse all <imputed_file> <empirical_file> -p 05_compare/ai2_vs_empirical/isec/

## Interpretation
# 0000.vcf = private to imputed
# 0001.vcf = private to wgrs
# 0002.vcf = imputed, shared
# 0003.vcf = wgrs, shared

# Save the shared output, and index
bcftools view 05_compare/ai2_vs_empirical/isec/0002.vcf -Ob -o 05_compare/ai2_vs_empirical/all_inds_imputed_shared.bcf

bcftools index 05_compare/ai2_vs_empirical/all_inds_imputed_shared.bcf 

bcftools view 05_compare/ai2_vs_empirical/isec/0003.vcf -Ob -o 05_compare/ai2_vs_empirical/all_inds_empirical_shared.bcf

bcftools index 05_compare/ai2_vs_empirical/all_inds_empirical_shared.bcf

# Clean space
rm -rf 05_compare/ai2_vs_empirical/isec
```

Compare the shared files:    
```
# Compare the genotypes between files, generate output files for comparison
./01_scripts/run_bcftools_stats.sh
# note: ...set user variables, including target folders and whether per-site should be calculated   

# Analyze the output, generate figures
01_scripts/assess_bcftools_stats.R
```


### 08. Calculate and evaluate impact of MAF ###
To inspect the effects of allele frequency or minor allele frequency, take the following steps:    
```
# Add all tags to the target imputed VCF file
bcftools +fill-tags 04_impute_all_loci/all_inds_wgrs_and_panel_biallele_only_fi3_imputed.vcf.gz -Ob -o 04_impute_all_loci/all_inds_wgrs_and_panel_biallele_only_fi3_imputed_with_tags.bcf
# ...here we will be interested in tags AF or MAF
# note: the MAF field will be the minor allele, not just the non-ref allele, so a completely homozygous alternate allele locus will have AF=1 and MAF=0

# If want to focus on AF for plotting alongside concordance from PSD file
bcftools query -f '%CHROM %POS %AF\n' file.bcf > output.txt 

# If want to focus on MAF
bcftools query -f '%CHROM %POS %MAF\n' 04_impute_all_loci/all_inds_wgrs_and_panel_biallele_only_fi3_imputed_with_tags.bcf > 04_impute_all_loci/all_inds_wgrs_and_panel_biallele_only_fi3_imputed_with_tags_perloc_MAF.txt

```

Plot with:    
`plot_psd_across_chr.R`    

If want to filter a BCF file and output only loci within a certain MAF range:    
```
bcftools view -i 'MAF>0 & MAF<0.05' eval_maf/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_offspring_only_rename_w_tags.bcf -Ob -o eval_maf/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_offspring_only_rename_w_tags_maf_0.001-0.5.bcf

```




