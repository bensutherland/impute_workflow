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
# Prepare to rename wgrs samples so they match between the files
bcftools query -l 11_impute_combine/parent_wgrs_only.vcf > 11_impute_combine/wgrs_parents_orgn_names.txt
# ...then manually annotate the file to separate the old_name and new_name with whitespace, single line per sample

# Rename in the VCF file
bcftools reheader --samples 11_impute_combine/wgrs_parents_orgn_names.txt -o 11_impute_combine/parent_wgrs_only_renamed.vcf 11_impute_combine/parent_wgrs_only.vcf

# Create sorted text file of new names
bcftools query -l 11_impute_combine/parent_wgrs_only_renamed.vcf | sort > 11_impute_combine/wgrs_parents_new_names_sorted.txt

# Sort in the VCF file
bcftools view -S 11_impute_combine/wgrs_parents_new_names_sorted.txt 11_impute_combine/parent_wgrs_only_renamed.vcf -o 11_impute_combine/parent_wgrs_only_renamed_sorted.vcf

# Compress and index 
bgzip 11_impute_combine/parent_wgrs_only_renamed_sorted.vcf

bcftools index 11_impute_combine/parent_wgrs_only_renamed_sorted.vcf.gz
```

Prepare the parent panel data to be combined:       
```
# Copy the parent panel data into the combined folder
cp -l 10_impute_input/parent_panel_roslin_rehead_hotspot_only.bcf 11_impute_combine/

# Prepare to rename panel samples so they match between the files
bcftools query -l 11_impute_combine/parent_panel_roslin_rehead_hotspot_only.bcf > 11_impute_combine/panel_parents_orgn_names.txt
# ...then manually annotate the file to separate the old_name and new_name with whitespace, single line per sample

# Rename in the VCF file
bcftools reheader --samples 11_impute_combine/panel_parents_orgn_names.txt 11_impute_combine/parent_panel_roslin_rehead_hotspot_only.bcf -o 11_impute_combine/parent_panel_roslin_rehead_hotspot_only_renamed.bcf

# Collect the new amp panel parent names and sort the names into a file
bcftools query -l 11_impute_combine/parent_panel_roslin_rehead_hotspot_only_renamed.bcf | sort > 11_impute_combine/panel_parents_new_names_sorted.txt

# Sort in the BCF file
bcftools view -S 11_impute_combine/panel_parents_new_names_sorted.txt 11_impute_combine/parent_panel_roslin_rehead_hotspot_only_renamed.bcf -o 11_impute_combine/parent_panel_roslin_rehead_hotspot_only_renamed_sorted.bcf

# index
bcftools index 11_impute_combine/parent_panel_roslin_rehead_hotspot_only_renamed_sorted.bcf

```

Combine the parent data with bcftools concat
```
bcftools concat --allow-overlaps 11_impute_combine/parent_wgrs_only_renamed_sorted.vcf.gz 11_impute_combine/parent_panel_roslin_rehead_hotspot_only_renamed_sorted.bcf -Ob -o 11_impute_combine/parent_wgrs_and_panel.bcf

# Index
bcftools index 11_impute_combine/parent_wgrs_and_panel.bcf

```

### 05. Merge parent wgrs and panel data with offspring panel data ###
a) Identify loci in the offspring panel file that overlap with the wgrs+panel parent file     
```
# Bring offspring data to folder
cp -l 10_impute_input/offspring_panel_roslin_rehead_hotspot_only.bcf* ./11_impute_combine/

# Create isec folder to capture output
mkdir 11_impute_combine/isec_combine_parents_and_offspring/

# Use isec to compare the files (no need for --collapse here, want only common shared REF alleles)
bcftools isec 11_impute_combine/parent_wgrs_and_panel.bcf 11_impute_combine/offspring_panel_roslin_rehead_hotspot_only.bcf -p 11_impute_combine/isec_combine_parents_and_offspring/

## Interpretation:    
# 0000.vcf = private to parents (wgrs+panel)
# 0001.vcf = private to offspring (panel)
# 0002.vcf = records from parents (wgrs+panel) shared in both
# 0003.vcf = records from offspring (panel) shared in both

# Save and rename 0003.vcf
cp 11_impute_combine/isec_combine_parents_and_offspring/0003.vcf 11_impute_combine/offspring_panel_roslin_rehead_hotspot_only_common_w_parents.vcf

# Compress and index
bgzip 11_impute_combine/offspring_panel_roslin_rehead_hotspot_only_common_w_parents.vcf
bcftools index 11_impute_combine/offspring_panel_roslin_rehead_hotspot_only_common_w_parents.vcf.gz
```

b) Combine the wgrs+panel parent data with the panel offspring data
```
bcftools merge 11_impute_combine/parent_wgrs_and_panel.bcf 11_impute_combine/offspring_panel_roslin_rehead_hotspot_only_common_w_parents.vcf.gz -Ob -o 11_impute_combine/all_inds_wgrs_and_panel.bcf

bcftools index 11_impute_combine/all_inds_wgrs_and_panel.bcf

# Copy the all-data file into the imputation folder, and index
cp -l 11_impute_combine/all_inds_wgrs_and_panel.bcf* 12_impute_impute/

# Remove multiallelic sites
bcftools view --max-alleles 2 ./12_impute_impute/all_inds_wgrs_and_panel.bcf -Ob -o 12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic.bcf    

```

#### Optional: Remove Mendelian incompatibility loci ####
Create a BCF file with parent and offspring panel-only loci:     
```
# Prepare an output folder for bcftools isec
mkdir 12_impute_impute/isec_keep_only_panel_loci/

# Index
bcftools index 12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic.bcf

# run isec to identify loci shared between all loci and panel-only offspring loci
bcftools isec ./12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic.bcf 11_impute_combine/offspring_panel_roslin_rehead_hotspot_only_common_w_parents.vcf.gz -p 12_impute_impute/isec_keep_only_panel_loci/

## Interpretation:    
# 0000.vcf = private to all_inds_wgrs_and_panel_no_multiallelic.bcf
# 0001.vcf = private to offspring panel
# 0002.vcf = records from all_inds shared in both
# 0003.vcf = records from offspring panel shared in both

# Save the records from all_inds shared in both
cp -l 12_impute_impute/isec_keep_only_panel_loci/0002.vcf 12_impute_impute/all_inds_panel_only.vcf

```

Use bcftools plugin mendelian to scan for Mendelian inconsistencies
```
# Use the pedigree file that has been annotated elsewhere in this pipeline, and format as needed for the plugin
awk '{ print $3 "," $2 "," $1 }' 12_impute_impute_no_novel/pedigree_annot.csv | grep -vE '^0' - > 12_impute_impute/pedigree.csv

# Use bcftools plugin Mendelian to annotate the number of Mendelian errors (MERR) in the BCF file and output
bcftools +mendelian 12_impute_impute/all_inds_panel_only.vcf -T 12_impute_impute/pedigree.csv --mode a -Ob -o 12_impute_impute/all_inds_panel_only_annot_MERR.vcf

# Observe the distribution of MERR
bcftools query -f '%CHROM %POS %MERR\n' 12_impute_impute/all_inds_panel_only_annot_MERR.vcf | sort -nk 3 | less

# Create a BCF file with the problematic loci 
bcftools view -i 'INFO/MERR >= 4' 12_impute_impute/all_inds_panel_only_annot_MERR.vcf -Ob -o 12_impute_impute/all_inds_panel_only_annot_MERR_problem_loci.bcf

# Index
bcftools index 12_impute_impute/all_inds_panel_only_annot_MERR_problem_loci.bcf
```

Remove the Mendelian inconsistencies from the wgrs+panel all individual file
```
# Prepare an output folder for bcftools isec
mkdir 12_impute_impute/isec_remove_MERR/

# run isec to identify loci private to the all loci data (dropping MERR)
bcftools isec ./12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic.bcf 11_impute_combine/offspring_panel_roslin_rehead_hotspot_only_common_w_parents.vcf.gz -p 12_impute_impute/isec_keep_only_panel_loci/

## Interpretation:    
# 0000.vcf = private to all_inds_wgrs_and_panel_no_multiallelic.bcf
# 0001.vcf = private to problem loci 
# 0002.vcf = records from all_inds shared in both
# 0003.vcf = records from problem loci shared in both

# Save the private records from all_inds
cp -l 12_impute_impute/isec_remove_MERR/0000.vcf 12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic_no_MERR.vcf

# Convert to BCF
bcftools view 12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic_no_MERR.vcf -Ob -o 12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic_no_MERR.bcf

# note: as clean-up, you may want to delete the isec folders, as they are large
```

#### Optional: Remove multi-mapper loci ####
Obtain the multi-mappers from Sutherland et al. 2024, Additional File S4. Save as two csv files, `bowtie_multimappers.csv` and `bwa_multimappers.csv` into `12_impute_impute`.    

Pull the loci from the latest BCF file that have marker names (i.e., they were from panel hotspots):    
```
bcftools view 12_impute_impute_no_novel_no_MERR/all_inds_wgrs_and_panel_no_multiallelic_no_MERR.bcf | grep -vE '^#' - | awk '$3!="." { print $0 }' - > 12_impute_impute/hotspot_loci_in_bcf.txt
# note: this file will be used to translate from the marker name to the chr and position info, which is the way the multimappers will be removed from the ai2 input file
```

Use Rscript to identify the chromosome and position names of the multimappers, read in the latest pre-impute ai2, and drop the multimappers from the pre-impute ai2 file:  
`01_scripts/impute_drop_multimappers_from_ai2.R`

This step will also require the ai2 input file that is produced from the BCF file (see next step), which will be subset using the specific markers and written back out. Then it joins the regular workflow by splitting the ai2 file into chromosomes and running ai2.    


#### Optional: Include some high density offspring to support the imputation ####
Copy in the all inds wgrs and panel BCF file, as well as the all offspring wgrs 10x data into `12_impute_impute`.     
Prepare files for combining:    
```
# Remove four offspring from each family from the imputation target file
bcftools view 12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic_no_MERR.bcf --samples ^ASY2_114_R1_1,ASY2_114_R1_2_ReAMP,ASY2_114_R1_3_ReAMP,ASY2_114_R1_4_ReAMP,ASY2_115_R1_1_ReAMP,ASY2_115_R1_2_ReAMP,ASY2_115_R1_3_ReAMP,ASY2_115_R1_4_ReAMP,ASY2_116_R1_1_ReAMP,ASY2_116_R1_2_ReAMP,ASY2_116_R2_1_ReAMP,ASY2_116_R2_2_ReAMP,ASY2_117_R1_1_ReAMP,ASY2_117_R1_2_ReAMP,ASY2_117_R1_3_ReAMP,ASY2_117_R1_4_ReAMP -Ob -o 12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic_no_MERR_rem_16_inds.bcf

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


### 08. Run GWAS on imputed data and plot results ###
Prepare GEMMA inputs:     
`01_scripts/imputed_ai2_to_gemma.R`      

Run GEMMA:    
```
cd 12_impute_impute
gemma -g gwas_geno.txt -p gwas_pheno.txt -gk -maf 0.05 -o gwas_all_fam
gemma -g gwas_geno.txt -p gwas_pheno.txt -k output/gwas_all_fam.cXX.txt -n 1 -c gwas_covar.txt  -maf 0.05 -lmm 4 -o gwas_all_fam_covar

```


Plot GEMMA outputs:    
`01_scripts/imputed_plot_gemma_results.R`    



### 09. Add grandparent data ###
Genotype grandparents with `wgrs_workflow`, then copy the filtered BCF to the present repository.     

```
# Index the two target files before running isec
bcftools index 12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic_no_MERR.bcf
bcftools index 12_impute_impute/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP10000_miss0.1.bcf

# prepare an output folder for bcftools isec
mkdir 12_impute_impute/combine_all_inds_and_grandparents/

# run isec 
bcftools isec 12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic_no_MERR.bcf 12_impute_impute/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP10000_miss0.1.bcf -p 12_impute_impute/combine_all_inds_and_grandparents/


## Interpretation:    
# 0000.vcf = private to all_inds_wgrs_and_panel_no_multiallelic_no_MERR.bcf
# 0001.vcf = private to grandparents (mpileup*) 
# 0002.vcf = records from all_inds shared in both
# 0003.vcf = records from grandparents shared in both

# Save out the target file to be combined
bcftools view 12_impute_impute/combine_all_inds_and_grandparents/0003.vcf -Ob -o 12_impute_impute/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP10000_miss0.1_compatible.bcf

# Index the output
bcftools index 12_impute_impute/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP10000_miss0.1_compatible.bcf

# Then can delete the isec folder to save space
```

Combine the all inds wgrs+panel with the grandparent RADseq     
```
bcftools merge 12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic_no_MERR.bcf 12_impute_impute/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP10000_miss0.1_compatible.bcf -Ob -o 12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic_no_MERR_w_grandparents.bcf
```

Next, go back up to the Imputation section and run.     


#### 06. Filter VCF and prepare for gemma analysis #### 
Use script `01_scripts/chr8_oshv1_amp_02_vcf_to_gemma.R`       


#### 07. Gemma analysis ####
Put output of the above script into `03_results`, then change directly into this folder to run the following commands.     
```
cd 03_results
gemma -g gwas_geno.txt -p gwas_pheno.txt -gk -maf 0.05 -o gwas_all_fam
gemma -g gwas_geno.txt -p gwas_pheno.txt -k output/gwas_all_fam.cXX.txt -n 1 -c gwas_covar.txt  -maf 0.05 -lmm 4 -o gwas_all_fam_covar
```

Then go to `chr8_oshv1_amp_03_gemma_results.R`.    

Optional: create a reduced dataset for testing:     
```
# Subset only a single chr for testing
bcftools view 11_impute_combine/all_inds_wgrs_and_panel.bcf --regions NC_047567.1 -Ob -o 12_impute_impute/all_inds_wgrs_and_panel_NC_047567_1.bcf

# Show number of panel loci
bcftools view 12_impute_impute/all_inds_wgrs_and_panel_NC_047567_1.bcf | grep -vE '^#' - | awk '{ print $3 }' - | sort | uniq -c | sort -nk1 | less

```        
