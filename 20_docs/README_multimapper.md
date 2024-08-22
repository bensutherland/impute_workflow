#### Optional: Remove multi-mapper loci ####
Obtain the multi-mappers from Sutherland et al. 2024, Additional File S4. Save as two csv files, `bowtie_multimappers.csv` and `bwa_multimappers.csv` into `12_impute_impute`.

Pull the loci from the latest BCF file that have marker names (i.e., they were from panel hotspots):
```
bcftools view 12_impute_impute_no_novel_no_MERR/all_inds_wgrs_and_panel_no_multiallelic_no_MERR.bcf | grep -vE '^#' - | awk '$3!="." { print $0 }' - > 12_impute_impute/hotspot_loci_in_bcf.txt
# note: this file will be used to translate from the marker name to the chr and position info, which is the way the multimappers will be removed from the ai2 input file
```

Use Rscript to identify the chromosome and position names of the multimappers, read in the latest pre-impute ai2, and drop the multimappers from the pre-impute ai2 file:
`01_scripts/impute_drop_multimappers_from_ai2.R`

This step will also require the ai2 input file that is produced from the BCF file (see next step), which will be subset using the specific markers and written back out. Then it joins the regular workflow by splitting the ai2 file into chromosomes and running ai2.
