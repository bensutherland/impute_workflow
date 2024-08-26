#!/bin/bash
# Run ai2 on individual chr prepared by 01_scripts/prep_geno_matrix_for_ai2.R
# B. Sutherland (2024-07-22)
# note: need to initiate ai2 in conda enviro first
# note: pedigree.csv is a space-delimited file

# Set user variables
INPUT_FOLDER="04_impute"
INPUT_PATTERN="ai2_input_"
PEDIGREE="pedigree.csv"

# Run ai2 iteratively
for file in $(ls -1 "$INPUT_FOLDER"/"$INPUT_PATTERN"*.txt)
do
    # Name of file
    echo "Imputation using $file"
    
    # basename
    name=$(basename "$file")

    # Impute 
    AlphaImpute2 -genotypes $INPUT_FOLDER/$name -pedigree $INPUT_FOLDER/$PEDIGREE -out $INPUT_FOLDER/${name%.txt} -maxthreads 12 -phase_output

done

#'2>&1' '>>' 10-log_files/"$TIMESTAMP"_01_cutadapt"${i%.fastq.gz}".log
