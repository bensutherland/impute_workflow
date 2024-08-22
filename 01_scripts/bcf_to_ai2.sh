#!/bin/bash
# Convert BCF file to ai2 genotype file 

# Set user variables
INPUT_FOLDER="04_impute"
INPUT_BCF="all_inds_wgrs_and_panel_biallele.bcf"

# Obtain vector of individual names from BCF file
echo "Preparing header for ai2 matrix"
bcftools query -l $INPUT_FOLDER/$INPUT_BCF > $INPUT_FOLDER/temp_ind_names.txt

# Add 'mname' at the top of vector, convert newlines to tabs, and remove last empty tab
echo "mname" | cat - $INPUT_FOLDER/temp_ind_names.txt |
        tr "\n" "\t" | 
        sed 's/\t$//g' > $INPUT_FOLDER/header.txt

# Add newline character to end of header file
echo "" >> $INPUT_FOLDER/header.txt

# Clean up space
rm $INPUT_FOLDER/temp_ind_names.txt

# BCFtools extract chr, position, and genotypes, then convert from genos to number alt alleles
echo "Preparing ai2 matrix from BCF file"
bcftools query -f '%CHROM %POS[\t%GT]\n' $INPUT_FOLDER/$INPUT_BCF | 
        sed 's/\b0\/0\b/0/g' | 
        sed 's/\b0\/1\b/1/g' | 
        sed 's/\b1\/1\b/2/g' | 
        sed 's/\.\/\./9/g' > $INPUT_FOLDER/${INPUT_BCF%.bcf}_ai2_no_header.txt

# Combine header with ai2 matrix
echo "Combining header with ai2 matrix"
cat $INPUT_FOLDER/header.txt $INPUT_FOLDER/${INPUT_BCF%.bcf}_ai2_no_header.txt > $INPUT_FOLDER/${INPUT_BCF%.bcf}_ai2.txt

# Remove temp file
rm $INPUT_FOLDER/${INPUT_BCF%.bcf}_ai2_no_header.txt
rm $INPUT_FOLDER/header.txt

