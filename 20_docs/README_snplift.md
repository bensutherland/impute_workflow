### 02. Convert to chromosome assembly coordinates ###
Move above the repo, then clone a snplift repo for each of offspring and parent data:
```
cd ..
git clone https://github.com/enormandeau/snplift.git snplift_offspring
git clone https://github.com/enormandeau/snplift.git snplift_parents

cp ./ms_cgig_chr8/10_impute_input/offspring_panel.vcf ./snplift_offspring/04_input_vcf/
cp ./ms_cgig_chr8/10_impute_input/parent_panel.vcf ./snplift_parents/04_input_vcf/
```

Within each snplift repository, edit `02_infos/snplift_config.sh`:
```
- full path to the original genome (bwa indexed)
- full path to the target genome (bwa indexed)
- relative path to the original VCF filename
- relative path to the new VCF filename
- CORRECT_ALLELES=1 to convert the ref all allele when alignments are reverse complemented
Note: setting the `CORRECT_ID` to 0 above prevents the ID column from being recalculated, so that your original IDs are carried through to the new VCF.
```

Run snplift for each:
```
cd snplift_offspring
time ./snplift 02_infos/snplift_config.sh
cp ./offspring_panel_roslin.vcf ../ms_cgig_chr8/10_impute_input/

cd ..

cd snplift_parents
time ./snplift 02_infos/snplift_config.sh
cp ./parent_panel_roslin.vcf ../ms_cgig_chr8/10_impute_input/

```

Return to `ms_cgig_chr8` repo, then further prepare panel data post-SNPlift:
```
# add headers
bcftools reheader 10_impute_input/offspring_panel_roslin.vcf --fai ~/genomes/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna.fai --output 10_impute_input/offspring_panel_roslin_rehead.vcf

bcftools reheader 10_impute_input/parent_panel_roslin.vcf --fai ~/genomes/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna.fai --output 10_impute_input/parent_panel_roslin_rehead.vcf

# compress
ls 10_impute_input/*_rehead.vcf | xargs -n 1 bgzip

# convert to BCF file
bcftools view 10_impute_input/offspring_panel_roslin_rehead.vcf.gz -Ob -o 10_impute_input/offspring_panel_roslin_rehead.bcf

bcftools view 10_impute_input/parent_panel_roslin_rehead.vcf.gz -Ob -o 10_impute_input/parent_panel_roslin_rehead.bcf

# index
ls 10_impute_input/*_rehead.bcf | xargs -n 1 bcftools index

```
