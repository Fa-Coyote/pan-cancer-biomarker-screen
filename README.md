# Pan-Cancer Candidate Biomarker Screening (Demonstration)

**Author:** Fatemeh Azizian Farsani

## What this is

A self-contained R script demonstrating a meta-analysis approach to
screening for genes that are consistently dysregulated across multiple
cancer types — the same general logic behind pan-cancer prognostic
biomarker discovery work I have contributed to as a co-author
(screening genes such as *FAM72B*, *CDK5R1*, *SERPINE1*, and *CISH*
across twelve cancer types).

**This script uses simulated, randomly generated expression data,
not real patient data or any third-party code.** It was written from
scratch to illustrate the method, since the original study's data and
code were produced collaboratively and are not mine to share.

## What it does

1. Simulates gene expression for Tumor vs. Normal samples across six
   "cancer types," with a handful of genes given a built-in
   differential signal and the rest left as noise.
2. Runs a per-cancer, per-gene differential expression test
   (two-sample t-test).
3. Combines the resulting p-values for each gene across all cancer
   types using **Fisher's method**, to find genes with a consistent
   signal across multiple cancer types rather than just one.
4. Applies FDR correction (Benjamini–Hochberg) and reports the
   top candidate biomarkers.

## How to run it

Requires only base R (no packages to install):

```bash
Rscript pan_cancer_biomarker_screen.R
```

## Background

This mirrors the analytical concept used in:

Shokrollahi A.\*, Mahdevar M.\*, Haji Ali Asgary Najafabadi A.H.,
**Azizian-Farsani F.**, Peymani M., Ghaedi K. *The expression of
FAM72B, MIR193BHG, CDK5R1, SERPINE1, and CISH as prognostic markers
in twelve prevalent cancers based on meta-analysis.* (Revised at
*Scientific Reports*)
