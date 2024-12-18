### Inspect empirical data for unexpected parental monomorphism ###
In addition to the Mendelian incompatibility section (see [here](20_docs/README_MI.md), additional information on genotype qualities can be obtained by observing instances where parents were fully genotyped but monomorphic, and the offspring were polymorphic.      

### 01. Prepare input data ###
Prepare the following files to be analyzed:    
```
# Parents-only renamed empirical VCF file
bcftools view 02_input_data/parent_hd/*_rename.bcf -Oz -o ./06_screen_loci/parents_only_rename.vcf.gz

# Offspring-only empirical VCF file
bcftools view 02_input_data/offspring_hd/*_rename.bcf -Oz -o 06_screen_loci/offspring_only_rename.vcf.gz

```

### 02. Analyze expected and unexpected monomorphism in parents per family ###
Use the following Rscript, which will have hard-coded family information:    
`01_scripts/shared_polymorphism_by_family.R`     
...output will be plots and data to screen. 


