# To do

## Taylor
- runtime SKAT vs our approach?
- Try MAF=0.001-0.01 with c=8,12, prevalence 5%
- Try prevalence<1%, like 0.1%
- Try craniofacial with MAF cutoff 3%
- Try simulations with 50kb
- table as opposed to qqplot

## Claudia
- simulation testing the proportional odds assumption: what happens if violated? -> added in text
- eigenvalue of M

# Associate Editor
Please note the request to address the proportional odds assumption, and the questions of prevalence and other assumptions in simulations in the impact of the method. I note that there was an attempt to marry case-control and linkage approaches by Hu et al in Pedigree VAAST (pVAAST) and am wondering if this should be referenced. 
- add reference to https://www.ncbi.nlm.nih.gov/pubmed/24837662

# Reviewer 2 (Comments for the Authors (Required)): 

The authors present a new approach to rare-variant association testing, called GAMuT, that leverages family history of disease to boost statistical power to detect associations. In particular they recommend the creation of an ordinal phenotype (0 for controls, 1 for sporadic cases, 2 for cases with a family history of disease) and derive appropriate statistics for testing aggregate rare-variant associations with this ordinal phenotype. Overall, the paper was well-written and logically structured. I have a few comments for the authors to consider. 

Comments 

(1) How important is the proportional odds assumption? Can the authors please provide some practical advice on how to determine whether results from GAMuT are still valid if the assumption is violated? i.e., how robust is GAMuT to violations of the proportional odds assumption? 
- relax proportional odds model assumption; we might want to do a simulation under no proportional odds and see how it performs

(2) Is run-time an issue? Some benchmarking against SKAT would be helpful. Some sequence-based datasets can be very large. 
- runtime benchmark vs SKAT: maybe they talk about this on the original gamut paper; did taylor simulations saving running time?

(3) For aggregate rare-variant associations, often the signal is due to just a few (4 or fewer) rare-variants. Does GAMuT help in those instances? I buy the rationale for an architecture where there is a gene with many very rare-risk variants (say 10 singletons). But examples of these types of architectures are pretty rare. More often we see a just a few intermediate frequency variants (0.001 - 0.01) with large-ish effect sizes. The simulations should be expanded to incorporate these more common architectures. 
- new simulations with MAF=0.001-0.01 with larger effect size: currently log10(MAF)*log(c)/4 with c=4,6
Try MAF=0.001 with c=8,12, prevalence 5%

(4) Though rare-variant association tests are often run for common diseases, they are also performed for rare diseases. It would be interesting to see the performance differences at an even smaller prevalence (< 1%). One would expect even more clustering of disease in families as the trait becomes more "Mendelian" in its architecture. 
- prevalence<1%, like 0.1% => might be good to run in AWS

Minor comment: 

(5) Why were the simulations conducted at a MAF cutoff of 3% whereas the data analysis was conducted at a MAF cutoff of 5%. It would look better if the thresholds were consistent or if the authors offered a rationale for this decision. 
- try craniofacial with MAF<3%


# Reviewer 3 (Comments for the Authors (Required)): 

Summary: 

The authors describe a natural extension of the standard case-control paradigm and a clever application of the existing software: GAMuT. They apply their proposed method to both simulated and real data. Their main result, that power to detect an association is increased when family history information is incorporated into the genetic analysis, appears to be sound... although I'm not quite sure why the variants need to be rare; some discussion on that point (if it's true) would be useful. 

Minor Comments: 

1. Researchers routinely exclude controls with a family history of disease because they (1) know/believe that the inclusion of such controls would reduce the power of a standard case-control approach, and (2) do not know of any software that would allow them to leverage information from such controls. (Furthermore, some researchers have both affected families and independent cases.) This time-honored practice is another line of evidence supporting the assertion that family-history matters, and that the author's method (esp. the version that contains 4 levels) could be valuable. Some additional discussion of these points might be useful. 
- add discussion point

2. In Introduction, the authors write: "..since more risk variants tend to segregate..." I know that we live in a variant-centric world, but variants do not segregate, FAMILIES segregrate! Thus, "...families with multiple affected individuals tend to segregrate more risk variants." is the correct statement. 
- change introduction

3. Using the cumulative-logit model to connect the ordinal problem to GAMuT was quite innovative. My only concern is that the constraint on the estimated multinomial probabilities (ie, summing to 1) should induce a co-linearity... thus, it would seem that the M matrix (and hence the I-M matrix) should both be dimension [Nx2]. A brief sentence or two on why [Nx3] is actually ok, or a re-write with [Nx2] would be appreciated. 
- check that there are only 2 eigenvalues on M, so that that takes care of the colinearity

4. The region size (eg 10kb) is a little small, given that the typical gene is between 50 and 100 kb in length. Some discussion on how ordinal GAMuT scales to larger regions would be useful. Also, wrt the same figure, Panel A, Panel B, Panel C, and Panel D, is much better than top, top middle, bottom middle, bottom.
- more simulations with 50kb

5. Albeit stylistic, some of us would rather see a type I error table with 10%, 5%, 1%, and 0.1%, since p-values > 5% are not so interesting. 
- add table to qqplots

6. Lastly, I know that the community has simply accepted that Li and Leal is the correct citation for burden, but in fact, the word burden cannot be found anywhere in the Li and Leal paper. As such, I think the authors should include a Madsen and Browning citation for the weighted sum statistic, which is much more in the spirit of burden. 
- add reference
