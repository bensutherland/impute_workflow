#!/bin/bash
# Run bcftools stats to compare two BCF files
# Assumes that you have a subfolder ('COMPARISON_FOLDER') within 'INPUT_FOLDER' that contains two indexed BCF files to be compared
# B. Sutherland (2024-09-17)

# Set user variables
INPUT_FOLDER="05_compare_all_loci"
COMPARISON_FOLDER="ai2_vs_empirical"
BCF1="all_inds_imputed_shared.bcf"
BCF2="all_inds_empirical_shared.bcf"
PER_SITE=true # set as true or false

## Compare the two BCF files
# Specify the samples to be run (set BCF2 as the file with fewest samples)
echo "Specifying samples to be compared in the two BCF files"
bcftools query -l $INPUT_FOLDER/$COMPARISON_FOLDER/$BCF2 >\
        $INPUT_FOLDER/$COMPARISON_FOLDER/inds_to_compare.txt

# Calculate stats on the two BCF files, for the specified individuals
if [ "$PER_SITE" = true ] ; then
  echo "Calculating stats on the two BCF files, including per-site discordance"
  bcftools stats --verbose -S $INPUT_FOLDER/$COMPARISON_FOLDER/inds_to_compare.txt \
          $INPUT_FOLDER/$COMPARISON_FOLDER/$BCF1 \
          $INPUT_FOLDER/$COMPARISON_FOLDER/$BCF2 >\
          $INPUT_FOLDER/$COMPARISON_FOLDER/stats_output.txt 

else  
  echo "Calculating stats on the two BCF files, without per-site discordance"
  bcftools stats -S $INPUT_FOLDER/$COMPARISON_FOLDER/inds_to_compare.txt \
          $INPUT_FOLDER/$COMPARISON_FOLDER/$BCF1 \
          $INPUT_FOLDER/$COMPARISON_FOLDER/$BCF2 >\
          $INPUT_FOLDER/$COMPARISON_FOLDER/stats_output.txt 
fi

## Collect required lines from output
#  Genotype Concordance Table
grep -E 'GCTs' $INPUT_FOLDER/$COMPARISON_FOLDER/stats_output.txt | 
        grep -v 'GCTs,' >\
         $INPUT_FOLDER/$COMPARISON_FOLDER/GCTs.txt

# Genotype concordance by sample for r-squared
grep -E 'GCsS' $INPUT_FOLDER/$COMPARISON_FOLDER/stats_output.txt | 
        grep -v 'GCsS,' >\
         $INPUT_FOLDER/$COMPARISON_FOLDER/GCsS.txt

# Per-site info (if collected)
if [ "$PER_SITE" = true ] ; then
  echo "Obtaining per-site information"
  grep -E 'PSD' $INPUT_FOLDER/$COMPARISON_FOLDER/stats_output.txt >\
  $INPUT_FOLDER/$COMPARISON_FOLDER/PSD.txt
fi
  
# Report finish and next step
echo "Stats calculated, sections subset and saved to output folder."
echo "Go to 01_scripts/assess_bcftools_stats.R"   

