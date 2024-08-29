#!/bin/bash
# Run Fimpute3 iteratively, must have used (...) #update
# B. Sutherland (2024-08-29)

# Set user variables
INPUT_FOLDER="04_impute/fimpute"
FIMPUTE="/home/greent/programs/FImpute3"

# Determine how many chr are present
NUM_CHR=$(ls -1 $INPUT_FOLDER/map_chr_* | sed 's/.*map_chr_//g' | sed 's/.txt//g' | wc -l)

# Change directory
cd $INPUT_FOLDER

# Prepare each control file
for i in $( seq 1 $NUM_CHR) 
do
   
   echo "Running FImpute, chr $i"
   
   $FIMPUTE "./fimpute_control_file_$i.txt" 

done


