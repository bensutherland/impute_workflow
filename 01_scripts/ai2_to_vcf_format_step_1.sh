#!/bin/bash
# Convert ai2 file to VCF file (step 1) 
# note: the INPUT_AI2 variable is the post-imputation file after it has been rejoined

# Set user variables
# AI2 settings
#INPUT_FOLDER="05_compare"
#INPUT_AI2="all_chr_combined.txt"

# FI3 settings
INPUT_FOLDER="04_impute/fimpute"
INPUT_AI2="fi3_loci_by_inds_all_imputed_chr.txt"

# Retain marker names, no header
awk -F"\t" '{ print $1 }' $INPUT_FOLDER/$INPUT_AI2 | 
        grep -vE '^mname' - > $INPUT_FOLDER/${INPUT_AI2%.txt}_mnames.txt

# Retain colnames 
head -n 1 $INPUT_FOLDER/$INPUT_AI2 > $INPUT_FOLDER/${INPUT_AI2%.txt}_header.txt

# Obtain only the genotypes (excluding first column and first row)
cut -f2- $INPUT_FOLDER/$INPUT_AI2 | tail -n +2 > $INPUT_FOLDER/${INPUT_AI2%.txt}_geno_only.txt

# Convert values back to VCF-style genotypes
sed 's/0/0\/0/g' $INPUT_FOLDER/${INPUT_AI2%.txt}_geno_only.txt | 
        sed 's/1/0\/1/g' - | 
        sed 's/2/1\/1/g' - > $INPUT_FOLDER/${INPUT_AI2%.txt}_geno_only_converted.txt

# Combine data back together (step 1)
paste -d "\t" $INPUT_FOLDER/${INPUT_AI2%.txt}_mnames.txt $INPUT_FOLDER/${INPUT_AI2%.txt}_geno_only_converted.txt > $INPUT_FOLDER/${INPUT_AI2%.txt}_converted_body.txt 

# Combine data back together (step 2)
cat $INPUT_FOLDER/${INPUT_AI2%.txt}_header.txt $INPUT_FOLDER/${INPUT_AI2%.txt}_converted_body.txt > $INPUT_FOLDER/${INPUT_AI2%.txt}_converted.txt

# Clean up
rm $INPUT_FOLDER/${INPUT_AI2%.txt}_mnames.txt
rm $INPUT_FOLDER/${INPUT_AI2%.txt}_header.txt
rm $INPUT_FOLDER/${INPUT_AI2%.txt}_geno_only.txt
rm $INPUT_FOLDER/${INPUT_AI2%.txt}_geno_only_converted.txt
rm $INPUT_FOLDER/${INPUT_AI2%.txt}_converted_body.txt

