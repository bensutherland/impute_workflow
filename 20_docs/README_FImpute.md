Instructions for preparing FImpute input from an ai2 prepared input file.     

To generate the map file, and the input materials for the genotypes file
R script: `01_scripts/ai2_to_fimpute3.R`     

Then go to terminal to finalize the genotypes file as per:    
```
# Convert 9 (ai2 missing) to 5 (FImpute missing)
sed 's/9/5/g' genotypes_no_indname.txt > genotypes_no_indname_9to5.txt

# Add header
echo 'Genotypes' | cat - genotypes_no_indname_9to5.txt > genotypes_no_indname_9to5_full.txt

# Combine all component parts
paste -d " " indnames.txt chip.txt genotypes_no_indname_9to5_full.txt  > genotype_file.txt

# Clean up workspace
rm genotypes_no_indname* indnames.txt chip.txt

```

Manually create a pedigree file in the following format:     
Animal Sire Dam Sex     


Then run FImpute3 and point to the control file.   Output will be in `output/genotypes_imp.txt`.    
