# impute_workflow
Imputation workflow for working with amplicon panel and wgrs data

### 00. Getting started ###
Clone the present repo, all commands will occur in the repo unless indicated.    

Requires the following inputs put in `02_input_data`, each in its own subfolder as labeled:      
- parent wgrs (20X) filtered genotypes BCF file from `wgrs_workflow` in `parent_hd`   
- parent amplicon panel individual vcf.gz or grouped bcf files in `parent_ld`   
- offspring amplicon panel individual vcf.gz or grouped bcf files `offspring_ld`   

Note: only include a single replicate per individual by the panel (pick the best replicate)   
Note: each type of amp panel data needed its own subfolder because they are named by the barcode and therefore will overlap with other datasets.     
Note: do not copy links of the VCF files, but rather full files.    


### 01. Prepare input data ### 
#### Separating and renaming the LD BCF file ####
If the input LD data is from amplitools _de novo_ genotyping of all parents and offspring together, it may be necessary to rename and separate the input BCF file. Put it in `02_input_data` root folder, and take the following steps (as an example):      
```
# Create a renaming file
bcftools query -l 02_input_data/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2.bcf > 02_input_data/sample_names_ld.csv
# ...then manually annotate the file by adding a space and giving the new name for each sample

# Rename samples in BCF file
bcftools reheader --samples 02_input_data/sample_names_ld.csv -o 02_input_data/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_rename.bcf 02_input_data/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2.bcf

# Separate parents and offspring LD data, here offspring noted by 'ASY2' in samplename
bcftools query -l 02_input_data/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_rename.bcf | grep -E '^ASY2' - > 02_input_data/offspring_samplelist.txt

bcftools query -l 02_input_data/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_rename.bcf | grep -vE '^ASY2' - > 02_input_data/parent_samplelist.txt

# Obtain offspring only
bcftools view -S 02_input_data/offspring_samplelist.txt 02_input_data/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_rename.bcf -Ob -o 02_input_data/offspring_ld/offspring_panel.bcf

bcftools index 02_input_data/offspring_ld/offspring_panel.bcf

# Obtain parents only
bcftools view -S 02_input_data/parent_samplelist.txt 02_input_data/mpileup_calls_noindel5_miss0.2_SNP_q0_avgDP10_biallele_minDP4_maxDP100000_miss0.2_rename.bcf -Ob -o 02_input_data/parent_ld/parent_panel.bcf

bcftools index 02_input_data/parent_ld/parent_panel.bcf
```


##### Combining offspring panel data #####
If the offspring data came as multiple VCF files, merge as follows:      
```
# decompress the files, then compress with bgzip
gunzip 10_impute_input/offspring_panel/*.gz
ls 10_impute_input/offspring_panel/*.vcf | xargs -n 1 bgzip

# index the files with bcftools
ls 10_impute_input/offspring_panel/*.vcf.gz | xargs -n 1 bcftools index

# create filelist for merging all VCF files
ls -1 10_impute_input/offspring_panel/*.vcf.gz > 10_impute_input/offspring_panel/sample_list.txt

# merge all VCF files
bcftools merge --file-list 10_impute_input/offspring_panel/sample_list.txt -Ov -o 10_impute_input/offspring_panel.vcf

```

##### Parent panel data #####
If the parent data came as multiple VCF files, merge as follows:       
```
# decompress the files, then compress with bgzip
gunzip 10_impute_input/parent_panel/*.gz
ls 10_impute_input/parent_panel/*.vcf | xargs -n 1 bgzip

# index the files with bcftools
ls 10_impute_input/parent_panel/*.vcf.gz | xargs -n 1 bcftools index

# create filelist for merging all VCF files
ls -1 10_impute_input/parent_panel/*.vcf.gz > 10_impute_input/parent_panel/sample_list.txt

# merge all VCF files
bcftools merge --file-list 10_impute_input/parent_panel/sample_list.txt -Ov -o 10_impute_input/parent_panel.vcf

```
note: these have novel variants also, which will need to be removed eventually (below)    


##### Parent wgrs data #####
If the parent and offspring wgrs data were all genotyped together, you will need to create a parent-only subset of the data, then index as follows:      

```
# Create a parent samplefile
bcftools query -l 10_impute_input/parent_wgrs/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1.bcf | grep -vE '^ASY2' - > 10_impute_input/parent_wgrs/parent_samplelist.txt

# Select only the parents from the dataset
bcftools view -S 10_impute_input/parent_wgrs/parent_samplelist.txt 10_impute_input/parent_wgrs/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1.bcf -Ob -o 10_impute_input/parent_wgrs/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_parents_only.bcf

# Index
bcftools index 10_impute_input/parent_wgrs/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_parents_only.bcf 
```


### 02. Exclude panel loci from parent wgrs data ###
Before merging the parent wgrs and panel loci, we need to remove all overlapping loci from the wgrs data. While doing this, it is possible to check for concordance in the overlapping loci.   

Separate loci into private or overlapping:     
```
# prepare an output folder for bcftools isec
mkdir 03_combine/isec_panel_and_wgrs

# run isec to identify overlapping or private loci. Include flag '--collapse all' to consider overlap regardless of alleles.    
bcftools isec --collapse all ./02_input_data/parent_hd/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_parents_only.bcf ./02_input_data/parent_ld/parent_panel.bcf -p 03_combine/isec_panel_and_wgrs/

## Interpretation:    
# 0000.vcf = private to parent wgrs
# 0001.vcf = private to parent panel
# 0002.vcf = records from parent wgrs shared in both
# 0003.vcf = records from parent panel shared in both

# Save the private records to parent wgrs
bcftools view 03_combine/isec_panel_and_wgrs/0000.vcf -Ob -o 03_combine/parent_wgrs_only.bcf

# Also save the overlapped loci for concordance evaluation
cp -l 03_combine/isec_panel_and_wgrs/0002.vcf 03_combine/parent_wgrs_shared_in_both.vcf
cp -l 03_combine/isec_panel_and_wgrs/0003.vcf 03_combine/parent_panel_shared_in_both.vcf

# Delete the isec folder to save space

```

### 03. Inspect concordance of shared loci in parents ###
- #TODO 


### 04. Concatenate parent panel loci into parent wgrs only data ###
Prepare the parent wgrs-only file to be combined:       
```
# Create a renaming file 
bcftools query -l 03_combine/parent_wgrs_only.bcf > 03_combine/parent_wgrs_samplenames.csv
# ...then manually annotate by adding a space and giving the new name for the sample 

# Rename in the BCF file
bcftools reheader --samples 03_combine/parent_wgrs_samplenames.csv -o 03_combine/parent_wgrs_only_renamed.bcf 03_combine/parent_wgrs_only.bcf

# Create sorted text file of new names
bcftools query -l 03_combine/parent_wgrs_only_renamed.bcf | sort > 03_combine/parent_wgrs_samplenames_sorted.csv

# Sort in the BCF file
bcftools view -S 03_combine/parent_wgrs_samplenames_sorted.csv 03_combine/parent_wgrs_only_renamed.bcf -Ob -o 03_combine/parent_wgrs_only_renamed_sorted.bcf

# Index 
bcftools index 11_impute_combine/parent_wgrs_only_renamed_sorted.vcf.gz
```

Prepare the parent panel data to be combined:       
```
# Copy the parent panel data into the combined folder
cp -l 02_input_data/parent_ld/parent_panel.bcf ./03_combine/

# No need to rename, but can follow the above instructions if needed

# Create sorted text file of names 
bcftools query -l 03_combine/parent_panel.bcf | sort > 03_combine/parent_panel_samplenames_sorted.csv

# Sort in the BCF file
bcftools view -S 03_combine/parent_panel_samplenames_sorted.csv 03_combine/parent_panel.bcf -Ob -o 03_combine/parent_panel_sorted.bcf

# index
bcftools index 03_combine/parent_panel_sorted.bcf
```

Combine the parent data with bcftools concat
```
bcftools concat --allow-overlaps 03_combine/parent_wgrs_only_renamed_sorted.bcf 03_combine/parent_panel_sorted.bcf -Ob -o 03_combine/parent_wgrs_and_panel.bcf

# Index
bcftools index 03_combine/parent_wgrs_and_panel.bcf
```


### 05. Merge parent wgrs and panel data with offspring panel data ###
If you need to identify which offspring panel data can be merged (i.e., if genotyping was not together), see #TODO README.     

Combine the wgrs+panel parent data with the panel offspring data
```
# Merge parent wgrs+panel with offspring panel
bcftools merge 03_combine/parent_wgrs_and_panel.bcf 02_input_data/offspring_ld/offspring_panel.bcf -Ob -o 03_combine/all_inds_wgrs_and_panel.bcf

# Remove multiallelic if present
bcftools view --max-alleles 2 03_combine/all_inds_wgrs_and_panel.bcf -Ob -o 03_combine/all_inds_wgrs_and_panel_biallele.bcf

# Copy to imputation folder
cp -l 03_combine/all_inds_wgrs_and_panel_biallele.bcf ./04_impute/

# Index 
bcftools index 04_impute/all_inds_wgrs_and_panel_biallele.bcf
```

- #TODO: add links to optional README approaches 


### 06. Imputation ###
The data is now all in a single BCF file and is ready for the imputation process.     

Prepare a pedigree file:      
```
bcftools query -l 12_impute_impute/<impute_target>.bcf > 12_impute_impute/pedigree.txt 

# Annotate the above file as follows:    
# <indiv> <sire> <dam>     
# where if there is no sire or dam, put 0
# save as space-delimited with the suffix `_annot.txt`.    
```

Format from BCF to AlphaImpute2:       
```
# Prepare ai2 matrix by building a header, and extracting info from the BCF, and converting to ai2 format
./01_scripts/bcf_to_ai2.sh
# produces: 12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic_ai2.txt

# the output will be in the following format:    
# mname \t ind1 \t ind2 \t (...)
# NC_047559.1 2945 \t 0 \t 0 \t 1 \t 9 (...)
# (...)
# where 0, 1, 2 are the number of alt alleles, and 9 is missing data  

# Split ai2 matrix into individual chr. This requires setting the input filename, and an identifiable chromosome string.
01_scripts/prep_geno_matrix_for_ai2.R   
# output will be in 12_impute_impute/ai2_input_<NC_047559.1>.txt, one file per chr 

```

Run imputation:     
```
# initialize the conda environment
conda activate ai2

# Run AlphaImpute2 on chromosome-separated datafiles
01_scripts/run_ai2.sh
# produces: 12_impute_impute/ai2_input_<NC_047559.1>.genotypes and *.haplotypes

# Transpose chromosome-separated imputed .genotypes files, and drop marker names on all matrices but the first to prepare for recombining the files back together
01_scripts/impute_rebuild_chr_lightweight.R
# produces: 13_impute_compare/*.genotypes_transposed_to_combine.txt

# Combine imputed, transposed, chr-separated files back together, then add marker names back in, based on the input ai2 file (before chromosome separation) 
#   note: before running, update the variable for your input original ai2 file
01_scripts/combine_transposed_ai2_output_and_mnames.sh
# produces: 13_impute_compare/all_chr_combined.txt

```

### 07. Evaluate imputation ###
Evaluate results by comparing the imputed data with the 10X 'empirical' data:     
```
# Obtain 10X bcf file 
cp -l ~/Documents/cgig/CHR8_wgrs/wgrs_workflow_offspring/05_genotyping/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1.bcf ./13_impute_compare/

# Use bash script to pull out genotypes into text file in ai2 format
# Edit the following script to point to the above bcf file, and run
01_scripts/bcf_to_ai2.sh

# Compare the imputed and empirical ai2 file genotypes by chromosome using: 
01_scripts/eval_impute_lightweight.R
# produces plots of average concordance between methods per individual by chromosome, and other outputs to screen
```

