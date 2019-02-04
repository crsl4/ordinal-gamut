README file for p_ACT.R and p_ACT_seq.R>
P-values Adjusted for Correlated Tests>

p_ACT.R is a post-processing program that adjusts groups of up to 
1000 correlated association tests for multiple testing.  Association 
tests may be based on K traits and M genotypes; not all K*M 
combinations need be considered.  Different genetic models involving 
the same genotype may also be incorporated by considering each model 
as a separate genotype.  

p_ACT_seq.R is similar to p_ACT.R except that while p_ACT.R adjusts 
only the best p-value from the group of tests, p_ACT_seq.R will 
sequentially adjust p-values in order of significance, stopping when 
p_ACT for a test is greater than alpha, the desired significance 
level, since sequential testing is no longer valid past this point.  
This method follows the procedure of Holm (1979). 

To run these programs, both R (http://www.r-project.org/) and the 
add-on package mvtnorm must be installed.  
( http://cran.r-project.org/src/contrib/Descriptions/mvtnorm.html ) 
If computing p_ACT for more than 100 tests, version 0.8-0 or later 
of 'mvtnorm' must be installed, as earlier versions only handle
<=100 tests properly.

Programs may be run interactively in UNIX and Windows environments by 
typing 'source(p_ACT.R)' from within an R session. In UNIX environments, 
batch mode is recommended; simply type 'R CMD BATCH p_ACT.R' at the prompt.

Output file: p_ACT.out (or p_ACT_seq.out for p_ACT_seq.R')

Input files: pvalues.txt, genotype.txt, traits.txt, covar.txt (optional), 
alpha.txt (optional, for p_ACT_seq.R only)
   If only a single trait is considered, the traits.txt file is not needed.
   If only a single marker is tested for association (with multiple traits), 
      the genotype.txt and covar.txt files are not needed.

Input file formats: space-delimited ASCII files.

pvalues.txt: 
   File should contain one row for each test, and either 3 or 4 columns.  
   Columns 1-2 should contain trait and SNP labels to distinguish the tests. 
   Column 3 should contain the p-value from each test.  
        No missing values allowed (so tests with missing p-values should be 
	   omitted from this file.)
   Column 4:
      If all tests are two-sided, no fourth column is needed.  
      If any one-sided tests were performed, Column 4 must contain:
         2 for two-sided tests 
	 Either -1 or 1 for one-sided tests, where sign indicates the 
	     direction of the test.  

genotype.txt: 
    One column for each genotype score. 
    One row for each individual plus a header row of genotype labels.  
    Missing values must be coded as NA.
    Not necessary if only one marker is considered.

covar.txt: 
    Not necessary if no covariates used.
    Not necessary if only one marker is considered.
    One column for each covariate. 
    One row for each individual plus header row of covariate names.
    Missing values coded as NA.  
    Must have same # of individuals as genotype.txt.

traits.txt: 
    Not necessary if only one trait is considered.
    One column for each trait. 
    One row for each individual plus header row of trait names.  
    If covariates are used, need trait values residualized on covariates.
    Missing values coded as NA.  
    Must have same # of individuals as genotype.txt.

alpha.txt
    Used only with p_ACT_seq.R (to determine when to stop sequential testing).
    Contains a single number indicating alpha level.
    If file not present, default value of .05 is used.

Example files for 7 tests on 6 individuals (K=3, M=4).  

pvalues.txt:
     bmi SNP1 .01 -1
     bmi SNP2_add .2  2
     bmi SNP2_dom .1  2
     bmi SNP2_rec .4  2
     dbp SNP1 .5  1
     dbp SNP2_add .1  2
     sbp SNP2_add .03  -1

 (Note that all traits are not tested on all genotypes; this is fine.)

genotype.txt: 
     SNP1 SNP2_add SNP2_dom SNP2_rec
     0 1 1 0
     0 2 1 1
     1 1 1 0
     2 1 1 0
     NA 0 0 0
     1 0 0 0

  (Note that SNP2 is run on three models, each having its own column in genotype.txt.)

covar.txt:
     c1 c2 c3
     7 2 3
     1 3 2
     1 NA 3
     1 5 3
     7 4 0
     0 3 4
  
  (If no covariates in model, this file is unnecessary.)

traits.txt:
     dbp sbp bmi
     5 6 7
     3 2 1
     3 5 5
     6 2 NA
     1 3 6
     4 3 1

alpha.txt: 
     .01
