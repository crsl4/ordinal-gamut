# Taylor project 

Meeting 8/27/2018

- Response: 0=control, 1=case without family history, 2=case with family history
- Families with parents, and different number of siblings (all probands have the same number of siblings)

Folder structure:
- Choo-documents: code and files from Choo, who was previous student (with Debashis) working on this. She described two approaches, we only used "logit-based approach)

- simulation-code: 
    - instructions.txt: instructions
    - folder have a name `pre` or `post` whether they have results or not
    - folder with `for-andrei` are run in AWS, others in HGCC
    - each `run` folder corresponds to one parameter setting, inside you have output files depending on whatever number of simulations you specified in `03` script
    - you specify the arguments in `02` script, and the number of simulations in `03`
    - `there is a spreadsheet with all the different folders, and the settings as columns:
        - sscheme: sampling scheme, not used
        - p adverse: probability that the causal marker is deleterious, not used

- results: results for power simulations (helper functions script with functions to get figures)
- inside `error-for-andrei-20180710` folder, we have the results for type I error
- `script.R` has the code to do the figures, and all output files are `Rmd`


- data-cleft: Rda file and bed/bim/fam, all input files come from Elizabeth
- `script.R` has the code to extract "good individuals" based on certain criteria: saved as `ids-to-keep.txt`. This also write the plink phenotype file with `FH` covariate, which is oue 0-1-2 response (explained above)
- `terminal script` is used to run plink to get a subset of the variants. you get `filtered` files and `vcf-filtered` (for bystro)
- the only covariates used are FH, PC1-5
- you upload the vcf files to bystro.io, which will identify to which genes each variant belongs to (annotation). The results are downloaded in `bystro` folder
- we read this folder (still in `script.R`), and write the `gene set.txt` file with a list of all the genes and which variants each one contains
- `script.R` (line 102): you can specify the minimum number of variants for a gene to be included in the analysis (here 4): `good_genes.set`: plink format file
- `set call script` creates one folder per gene, and creates the raw files with plink (which are the 0-1-2 matrices of minor allele counts, which we need for R)
- the function `do_gene_test` does the multinomial regression per gene. Note that it should use the "good genes" as input, not the `rev_good_genes`. Taylor had to create the `rev_good_genes` because of code that stopped
- This function will produce `gene-test-results.txt` with the pvalues of three tests: burden, skat and gamut (ours)
- then we added location info to his file: `res_with_loc.Rda`

- folder `results_cleft`: Rmd and html of results

# Other resources
  - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4800053/
  - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4833215/

# Rmd manuscript
- Created the `ordinal-skat.Rmd` file
- Had to `brew install pandoc-citeproc` to use bib file, and `brew install pandoc`.
- Had to make sure all these packages were installed in R:
```r
xfun::session_info(c(
  'blogdown', 'bookdown', 'knitr', 'rmarkdown', 'htmltools',
  'reticulate', 'rticles', 'flexdashboard', 'learnr', 'shiny',
  'revealjs', 'pkgdown', 'tinytex', 'xaringan', 'tufte'
), dependencies = FALSE)
```
with the versions greater than or equal to [here](https://bookdown.org/yihui/rmarkdown/software-info.html)

- Power plots from `ordinal_results_parents_addtl_variant_stats_20180605.html`:
```shell
cd Dropbox/Documents/gwas/projects/ordinal-skat/scripts
scp csolislemus@hgcc.genetics.emory.edu:/home/csolislemus/ordinal_family_history_files_for_claudia/simulation_code/results/ordinal_results_parents_addtl_variant_stats_20180605.* .
```

Then, I copied as `ordinal_results_parents_addtl_variant_stats_20180605-modified.Rmd`, which I will modify for the specific power plots that we need: columns (lambda=2,4,8), rows (prevalence=0.01,0.05), this matrix for one effect size, and results only for parents as family. We also use this Rmd for the causal variants plots.
This Rmd file has the code to create specific CSV files that will later be used inside `ordinal-gamut.Rmd` for the plots.


# Running bystro
Hi Claudia,

No worries! Hope you are having a nice rainy-but-not-really-snowy day :) Here’s an outline for running Bystro:

1. When you run “terminal_script”, the last line generates “vcf_filtered.vcf”
2. Go to https://bystro.io/
3. Click start and log in
4. Select Hg19 assembly
5. Upload the vcf file
6. Once the upload is completed, go to the results tab
7. Where it says “search this file”, enter “-intergenic”
8. Once filter is done, click the little blue floppy disk icon to save the query results
9. Click on results tab again
10. Select the annotation that used the -intergenic query
11. Click the icon with 3 blue lines and download annotation
12. Open the downloaded file ending in “annotation.tsv.gz” to generate the “annotation.tsv” file
13. You can then continue running script.R, which will read in and format that file

Let me know if you have any problems!

Best,
Taylor