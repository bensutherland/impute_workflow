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
- Screen for Mendelian incompatibilities in panel data [here](20_docs/README_MI.mdi)     


### 01. Prepare input data ### 
#### wgrs data (HD data) ####
##### Separate and rename parents and offspring #####   
Separate and rename parents:      
```
# Create a parent samplefile
bcftools query -l 02_input_data/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1.bcf | grep -vE '^ASY2' - > 02_input_data/parent_hd_samplelist.txt
# ...note: if there are any other parents you want to remove from the dataset, delete them from this select list here.     

# Select only the parents from the dataset
bcftools view -S 02_input_data/parent_hd_samplelist.txt 02_input_data/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1.bcf -Ob -o 02_input_data/parent_hd/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_parents_only.bcf

# Rename 
# ...manually add desired samplenames to 02_input_data/parent_hd_samplelist.txt in space-separated format
bcftools reheader --samples 02_input_data/parent_hd_samplelist.txt -o 02_input_data/parent_hd/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_parents_only_rename.bcf 02_input_data/parent_hd/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_parents_only.bcf

# Index
bcftools index 02_input_data/parent_hd/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_parents_only_rename.bcf

# Clean space
rm 02_input_data/parent_hd/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_parents_only.bcf
```

Separate and rename offspring:      
```
# Create an offspring samplefile
bcftools query -l 02_input_data/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1.bcf | grep -E '^ASY2' - > 02_input_data/offspring_hd_samplelist.txt
# ...note: if there are any other offspring you want to remove from the dataset, delete them from this select list here

# Select only the offspring from the dataset
bcftools view -S 02_input_data/offspring_hd_samplelist.txt 02_input_data/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1.bcf -Ob -o 02_input_data/offspring_hd/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_offspring_only.bcf

# Rename
# ...manually add desired samplenames to the offspring samplefile in space-sep format
awk -F"-" '{ print $0 " " $1"-"$2"-"$3"-"$4 }' 02_input_data/offspring_hd_samplelist.txt > 02_input_data/offspring_hd_samplelist_for_rename.txt

bcftools reheader --samples 02_input_data/offspring_hd_samplelist_for_rename.txt -o 02_input_data/offspring_hd/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_offspring_only_rename.bcf 02_input_data/offspring_hd/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_offspring_only.bcf

# Index
bcftools index 02_input_data/offspring_hd/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_offspring_only_rename.bcf

# Clean space
rm 02_input_data/offspring_hd/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_offspring_only.bcf
```

#### panel data (LD data) ####
##### Separate and rename parents and offspring #####   
Separate and rename parents:      
```
# Create a parent samplefile (in this example, identified by string 'rawlib')
bcftools query -l 02_input_data/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2.bcf | grep 'rawlib' - > 02_input_data/parent_ld_samplelist.txt
# ...note: if there are any parents you want to remove, delete from this samplelist

# Select only the parents from the dataset
bcftools view -S 02_input_data/parent_ld_samplelist.txt 02_input_data/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2.bcf -Ob -o 02_input_data/parent_ld/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_parents_only.bcf

# Rename
# ...manually add desired samplenames to 02_input_data/parent_ld_samplelist.txt in space-separated format
bcftools reheader --samples 02_input_data/parent_ld_samplelist.txt -o 02_input_data/parent_ld/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_parents_only_rename.bcf 02_input_data/parent_ld/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_parents_only.bcf

# Index
bcftools index 02_input_data/parent_ld/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_parents_only_rename.bcf
```

Separate and rename offspring:    
```
# Create an offspring samplefile (in this example, identified by absence of 'rawlib')
bcftools query -l 02_input_data/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2.bcf | grep -vE 'rawlib' - > 02_input_data/offspring_ld_samplelist.txt
# ...note: if there are any offspring you want to remove, delete from this samplelist

# Select only the offspring from the dataset
bcftools view -S 02_input_data/offspring_ld_samplelist.txt 02_input_data/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2.bcf -Ob -o 02_input_data/offspring_ld/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_offspring_only.bcf

# Rename
# ...manually create 02_input_data/offspring_ld_samplelist_for_rename.txt
bcftools reheader --samples 02_input_data/offspring_ld_samplelist_for_rename.txt -o 02_input_data/offspring_ld/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_offspring_only_rename.bcf 02_input_data/offspring_ld/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_offspring_only.bcf

# Index
bcftools index 02_input_data/offspring_ld/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_offspring_only_rename.bcf
```


### 02. Exclude panel loci from parent wgrs data ###
Before merging parent wgrs and panel data, all overlapping loci between the two will be dropped from the wgrs dataset. Note: this will also give substrate for direct comparison of overlapping loci (see next section).     

Separate loci into private or overlapping:     
```
# prepare an output folder for bcftools isec
mkdir 03_combine/isec_wgrs_and_panel

# run isec to identify overlapping or private loci. Include flag '--collapse all' to consider overlap regardless of alleles.    
bcftools isec --collapse all 02_input_data/parent_hd/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_parents_only_rename.bcf 02_input_data/parent_ld/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_parents_only_rename.bcf -p 03_combine/isec_wgrs_and_panel/

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
cp -l 02_input_data/parent_ld/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_parents_only_rename.bcf 03_combine/

# Create sorted text file of names 
bcftools query -l 03_combine/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_parents_only_rename.bcf | sort > 03_combine/parent_panel_samplenames_sorted.txt

# Sort in the BCF file
bcftools view -S 03_combine/parent_panel_samplenames_sorted.txt 03_combine/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_parents_only_rename.bcf -Ob -o 03_combine/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_parents_only_rename_sorted.bcf

# index
bcftools index 03_combine/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_parents_only_rename_sorted.bcf
```

Combine the parent data with bcftools concat
```
bcftools concat --allow-overlaps 03_combine/parent_wgrs_only_sorted.bcf 03_combine/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_parents_only_rename_sorted.bcf -Ob -o 03_combine/parent_wgrs_and_panel.bcf

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

Prepare a pedigree file:      
```
bcftools query -l 04_impute/all_inds_wgrs_and_panel_biallele.bcf > 04_impute/pedigree.csv

# Annotate the above file as follows:    
# <indiv> <sire> <dam>     
# where if there is no sire or dam, put 0
# save as space-delimited    
```

Imputation with AlphaImpute2:       
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

Imputation with FImpute3:    
```
# Note: the version of FImpute3 used here can only handle fewer than 2 M SNPs per run, so it is necessary to separate the input file into chromosomes before running, then run FImpute3 per chr (as with ai2).     

# Note: the first step expects that the BCF file has already been converted to ai2 format 

# Convert from ai2 to FImpute3 format
./01_scripts/ai2_to_fimpute3.R
# ...this creates a pedigree file, a map file per chr, and parts needed for genotypes file per chr

# Edit and combine components needed for per-chr genotypes files
./01_scripts/prep_fi3_indiv_chr.sh

# Prepare a control file and output run directory per chr 
./01_scripts/prep_fi3_control_files_and_output_dir.sh

# Run FImpute3 iteratively
./01_scripts/run_fi3_iteratively.sh

# Convert the FImpute3 output back to alphaimpute2 format for evaluation
./01_scripts/fimpute3_to_ai2.R

```



### 07. Evaluate imputation ###
Evaluate results by comparing the imputed data with the 10X 'empirical' data:     
```
# Obtain 10X bcf file 
cp -l 02_input_data/offspring_hd/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_offspring_only_rename.bcf ./05_compare/

# Convert from BCF to ai2 after updating the path and filename
01_scripts/bcf_to_ai2.sh

# Compare the imputed and empirical ai2 file genotypes by chromosome using: 
01_scripts/eval_impute_lightweight.R
# produces plots of average concordance between methods per individual by chromosome, and other outputs to screen
```


##### New approach: #####    
```
# Identify which regions were output by ai2
awk '{ print $1 "\t" $2  }' 05_compare/all_chr_combined.txt | grep -vE '^mname' - > 05_compare/ai2_imputed_regions.txt

# Prepare a VCF file from the pre-imputation BCF file, with only the selected regions
bcftools view --regions-file 05_compare/ai2_imputed_regions.txt 04_impute/all_inds_wgrs_and_panel_biallele.bcf -Ov -o 04_impute/all_inds_wgrs_and_panel_biallele_only_ai2_imputed_regions.vcf

# Format the ai2 output file to prepare it to be used to create a VCF file
./01_scripts/ai2_to_vcf_format_step_1.sh

# Use the following Rscript to read in the VCF and fill in the data from the imputed file
01_scripts/ai2_to_VCF.R

```

Compare the imputed to the empirical:      



### 08. Inspect concordance of shared loci in parents and offspring ###
Evaluate the concordance of the shared loci in the panel and wgrs data for both the parents and offspring. First, prepare the data:     
```
# Make a subfolder to keep things tidy
mkdir 05_compare/panel_vs_wgrs

# Copy in all datasets (note: filenames must be different)   
cp -l 02_input_data/parent_hd/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_parents_only_rename.bcf* 05_compare/panel_vs_wgrs/

cp -l 02_input_data/parent_ld/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_parents_only_rename.bcf* ./05_compare/panel_vs_wgrs/

cp -l 02_input_data/offspring_hd/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_offspring_only_rename.bcf* 05_compare/panel_vs_wgrs/

cp -l 02_input_data/offspring_ld/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_offspring_only_rename.bcf* 05_compare/panel_vs_wgrs/

# Merge the wgrs data back together
bcftools merge 05_compare/panel_vs_wgrs/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_parents_only_rename.bcf 05_compare/panel_vs_wgrs/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_offspring_only_rename.bcf -Ob -o 05_compare/panel_vs_wgrs/parents_and_offspring_wgrs_renamed.bcf  

bcftools index 05_compare/panel_vs_wgrs/parents_and_offspring_wgrs_renamed.bcf

# Merge the panel data back together
bcftools merge 05_compare/panel_vs_wgrs/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_parents_only_rename.bcf 05_compare/panel_vs_wgrs/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_offspring_only_rename.bcf -Ob -o 05_compare/panel_vs_wgrs/parents_and_offspring_panel_renamed.bcf

bcftools index 05_compare/panel_vs_wgrs/parents_and_offspring_panel_renamed.bcf

```

Prepare files that only contain common loci between the platforms:     
```
mkdir 05_compare/panel_vs_wgrs/isec/

# run isec to identify overlapping or private loci. Include flag '--collapse all' to consider overlap regardless of alleles.    
bcftools isec --collapse all 05_compare/panel_vs_wgrs/parents_and_offspring_wgrs_renamed.bcf 05_compare/panel_vs_wgrs/parents_and_offspring_panel_renamed.bcf -p 05_compare/panel_vs_wgrs/isec/

## Interpretation:    
# 0000.vcf = wgrs, private
# 0001.vcf = panel, private
# 0002.vcf = wgrs, shared
# 0003.vcf = panel, shared

# Save the shared output, this will be used for comparisons
bcftools view 05_compare/panel_vs_wgrs/isec/0002.vcf -Ob -o 05_compare/panel_vs_wgrs/all_inds_wgrs_shared.bcf  

bcftools index 05_compare/panel_vs_wgrs/all_inds_wgrs_shared.bcf

bcftools view 05_compare/panel_vs_wgrs/isec/0003.vcf -Ob -o 05_compare/panel_vs_wgrs/all_inds_panel_shared.bcf

bcftools index 05_compare/panel_vs_wgrs/all_inds_panel_shared.bcf

# Clean the space
rm -rf 05_compare/panel_vs_wgrs/isec
```

Make the actual comparisons
```
# Specify the samples that you want to run (must be present in both files)
bcftools query -l 05_compare/panel_vs_wgrs/all_inds_wgrs_shared.bcf > 05_compare/panel_vs_wgrs/inds_to_compare.txt
# ...note: do this on the wgrs data (the datatype with the fewer samples)

bcftools stats -S 05_compare/panel_vs_wgrs/inds_to_compare.txt  05_compare/panel_vs_wgrs/all_inds_wgrs_shared.bcf 05_compare/panel_vs_wgrs/all_inds_panel_shared.bcf > 05_compare/panel_vs_wgrs/all_inds_wgrs_panel_comp.txt

grep -E 'GCTs' 05_compare/panel_vs_wgrs/all_inds_wgrs_panel_comp.txt  | grep -v 'GCTs,' > 05_compare/panel_vs_wgrs/all_inds_wgrs_panel_comp_GCTs.txt

grep -E 'GCsS' 05_compare/panel_vs_wgrs/all_inds_wgrs_panel_comp.txt  | grep -v 'GCsS,' - > 05_compare/panel_vs_wgrs/all_inds_wgrs_panel_comp_GCsS.txt
```

Use the following script to analyze the outputs and generate figures:    
`01_scripts/assess_bcftools_stats.R`    












 
