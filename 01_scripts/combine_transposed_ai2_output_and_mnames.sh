#!/bin/bash
# Combine transposed output and add mname column back into data
# B. Sutherland (2024-07-23)

# Set user variables
INPUT_FOLDER="05_compare/"
INPUT_PATTERN=".genotypes_transposed_to_combine.txt"
ORGN_AI2_FILE="04_impute/all_inds_wgrs_and_panel_biallele_ai2.txt"

# Combine all imputed datafiles 
cat "$INPUT_FOLDER"/*"$INPUT_PATTERN" > "$INPUT_FOLDER"/all_chr_combined_temp.txt

# Collect mname column
awk -F"\t" '{ print $1 }' $ORGN_AI2_FILE |
        grep -E 'mname|NC_047' - > "$INPUT_FOLDER"/mnames_temp.txt

# Add mname
paste $INPUT_FOLDER/mnames_temp.txt "$INPUT_FOLDER"/all_chr_combined_temp.txt > "$INPUT_FOLDER"/all_chr_combined.txt

# Clean up
rm "$INPUT_FOLDER"/all_chr_combined_temp.txt
rm "$INPUT_FOLDER"/mnames_temp.txt

