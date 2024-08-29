#!/bin/bash
# Prepare FImpute3 file, must already have run 01_scripts/ai2_to_fimpute3.R
# B. Sutherland (2024-08-29)

# Set user variables
INPUT_FOLDER="04_impute/fimpute"

# Determine how many chr are present
NUM_CHR=$(ls -1 $INPUT_FOLDER/map_chr_* | sed 's/.*map_chr_//g' | sed 's/.txt//g' | wc -l)

# Prepare each chromosome
for i in $( seq 1 $NUM_CHR) 
do
   
   echo "Working on chr $i"

   # Convert 9s (missing val for ai2) to 5s (missing val for fi3)
   sed 's/9/5/g' $INPUT_FOLDER/"genotypes_no_indname_chr_$i.txt" > $INPUT_FOLDER/"genotypes_no_indname_9to5_chr_$i.txt"   
   
   # Add header
   echo 'Genotypes' | cat - $INPUT_FOLDER/"genotypes_no_indname_9to5_chr_$i.txt" > $INPUT_FOLDER/"genotypes_no_indname_9to5_full_chr_$i.txt"
  
   # Combine all component parts 
   paste -d " " $INPUT_FOLDER/indnames.txt $INPUT_FOLDER/chip.txt $INPUT_FOLDER/"genotypes_no_indname_9to5_full_chr_$i.txt" > $INPUT_FOLDER/"genotype_file_chr_$i.txt" 

done


  
 # Clean up workspace
 rm $INPUT_FOLDER/genotypes_no_indname_9to5* 
 # 
 # # Reporting
 # echo "Results are now available in $INPUT_FOLDER/genotype_file.txt"
 # 
