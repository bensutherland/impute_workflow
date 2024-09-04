#!/bin/bash
# Convert ai2 file to VCF file (step 1) 

# Set user variables
INPUT_FOLDER="05_compare"
INPUT_AI2="all_chr_combined.txt"

# Retain marker names
awk -F"\t" '{ print $1 }' 05_compare/all_chr_combined.txt > 05_compare/all_chr_combined_mnames.txt

# Retain colnames 
#TODO#

# Obtain only the genotypes (excluding first column
cut -f2- 05_compare/all_chr_combined.txt > 05_compare/all_chr_combined_geno_only.txt

# Convert values back to VCF-style genotypes
sed 's/0/0\/0/g'  05_compare/all_chr_combined_geno_only.txt | sed 's/1/0\/1/g' - | sed 's/2/1\/1/g' - > 05_compare/all_chr_combined_geno_only_converted.txt

# Recombine
paste -d "\t" 05_compare/all_chr_combined_mnames.txt 05_compare/all_chr_combined_geno_only_converted.txt

