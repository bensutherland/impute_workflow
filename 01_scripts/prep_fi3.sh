#!/bin/bash
# Prepare FImpute3 file, must already have run 01_scripts/ai2_to_fimpute3.R
# note: ** no longer used **, replaced by 01_scripts/prep_fi3_indiv_chr.sh
# B. Sutherland (2024-08-27)

# Set user variables
INPUT_FOLDER="04_impute/fimpute/"

# Convert 9s (missing val for ai2) to 5s (missing val for fi3)
sed 's/9/5/g' $INPUT_FOLDER/genotypes_no_indname.txt > $INPUT_FOLDER/genotypes_no_indname_9to5.txt

# Add header to the genotypes section
echo 'Genotypes' | cat - $INPUT_FOLDER/genotypes_no_indname_9to5.txt > $INPUT_FOLDER/genotypes_no_indname_9to5_full.txt

# Combine all component parts 
paste -d " " $INPUT_FOLDER/indnames.txt $INPUT_FOLDER/chip.txt $INPUT_FOLDER/genotypes_no_indname_9to5_full.txt  > $INPUT_FOLDER/genotype_file.txt

# Clean up workspace
rm $INPUT_FOLDER/genotypes_no_indname* $INPUT_FOLDER/indnames.txt $INPUT_FOLDER/chip.txt

# Reporting
echo "Results are now available in $INPUT_FOLDER/genotype_file.txt"

