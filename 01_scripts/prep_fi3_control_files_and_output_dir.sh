#!/bin/bash
# Prepare FImpute3 control files, must already have run 01_scripts/ai2_to_fimpute3_.R
# B. Sutherland (2024-08-29)

# Set user variables
INPUT_FOLDER="04_impute/fimpute"
CONTROL_FILE_TEMPLATE="00_archive/fimpute3_control_file_template.txt"

# Determine how many chr are present
NUM_CHR=$(ls -1 $INPUT_FOLDER/map_chr_* | sed 's/.*map_chr_//g' | sed 's/.txt//g' | wc -l)

# Prepare each control file
for i in $( seq 1 $NUM_CHR) 
do
   
   echo "Working on control file $i"
   
   # Make control file per chr
   sed "s/genotype_file.txt/genotype_file_chr_$i.txt/g" $CONTROL_FILE_TEMPLATE  | 
       sed "s/map.txt/map_chr_$i.txt/g" | 
       sed "s/\"output\"/\"output_folder_$i\"/g" > $INPUT_FOLDER/fimpute_control_file_$i.txt 

   # Make output folders
   mkdir $INPUT_FOLDER/output_folder_$i

done

# Run Fimpute3 using the next script:  

