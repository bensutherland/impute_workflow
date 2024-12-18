### 09. Add grandparent data ###
Genotype grandparents with `wgrs_workflow`, then copy the filtered BCF to the present repository. Also copy in the renamed offspring and parents all empirical data from the standard workflow.     


```
# Index the two target files before running isec

# prepare an output folder for bcftools isec
mkdir 02_input_data/isec_wgrs_v_rad

# run isec
bcftools isec --collapse none 02_input_data/parents_and_offspring_wgrs.bcf 02_input_data/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP10000_miss0.1.bcf -p 02_input_data/isec_wgrs_v_rad/

## Interpretation:
# 0000.vcf = private to wgrs parents and offspring
# 0001.vcf = private to RAD grandparents
# 0002.vcf = records from wgrs shared in both
# 0003.vcf = records from RAD shared in both

# Save out the target file to be combined
bcftools view 02_input_data/isec_wgrs_v_rad/0002.vcf -Ob -o 02_input_data/parents_and_offspring_wgrs_shared_loci_w_grandparents.bcf

bcftools view 02_input_data/isec_wgrs_v_rad/0003.vcf -Ob -o 02_input_data/grandparents_shared_loci_w_parents_and_offspring.bcf

# Index the output

# Merge
bcftools merge 02_input_data/grandparents_shared_loci_w_parents_and_offspring.bcf 02_input_data/parents_and_offspring_wgrs_shared_loci_w_grandparents.bcf -Ob -o 02_input_data/grandparents_parents_and_offspring.bcf

```

Next, prepare for making an in silico panel, but first separate grandparent+parent from offspring:    
```
## Create sample lists for all three generations
# Offspring
bcftools query -l 02_input_data/grandparents_parents_and_offspring.bcf | grep -E '^ASY2' - > 02_input_data/offspring_samplelist.txt

# Parents
bcftools query -l 02_input_data/grandparents_parents_and_offspring.bcf | grep -vE '^SRR|^ASY2' - > 02_input_data/parents_samplelist.txt

# Grandparents
bcftools query -l 02_input_data/grandparents_parents_and_offspring.bcf | grep -E '^SRR' - > 02_input_data/grandparents_samplelist.txt

## Select only the offspring from the file 
bcftools view -S 02_input_data/offspring_samplelist.txt 02_input_data/grandparents_parents_and_offspring.bcf -Ob -o 02_input_data/offspring_all_loci.bcf

## Select the parents and grandparents from the file
cat 02_input_data/grandparents_samplelist.txt 02_input_data/parents_samplelist.txt > 02_input_data/grandparents_and_parents_samplelist.txt

bcftools view -S 02_input_data/grandparents_and_parents_samplelist.txt 02_input_data/grandparents_parents_and_offspring.bcf -Ob -o 02_input_data/grandparents_and_parents_all_loci.bcf

```

Go to README for in silico panel to create from the offspring BCF.    

Once that has been made, copy it into offspring LD folder:     
`cp 10_in_silico_panel/subsample_sort.bcf* 02_input_data/offspring_ld/`      

Combine the in silico panel data and the full data:     
(note: this section is similar to '04. Merge parent wgrs and panel data with offspring panel data' from the main README).       
```
bcftools merge 02_input_data/grandparents_and_parents_all_loci.bcf 02_input_data/offspring_ld/subsample_sort.bcf -Ob -o 03_combine/all_inds_wgrs_and_panel.bcf

# Copy to imputation folder
cp -l 03_combine/all_inds_wgrs_and_panel.bcf 04_impute/
# Index
bcftools index 04_impute/all_inds_wgrs_and_panel.bcf 
```

Go to the main README, and follow the 05. Imputation section.    




