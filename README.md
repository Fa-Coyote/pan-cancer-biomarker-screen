# Pan-Cancer Candidate Biomarker Screening (Demonstration)

**Author:** Fatemeh Azizian Farsani

## What this is

A self-contained R script demonstrating a meta-analysis approach to
screening for genes that are consistently dysregulated across multiple
cancer types. This mirrors the general logic behind a pan-cancer
prognostic biomarker discovery study I co-authored, currently under
journal review and not yet public — so this demo intentionally uses
simulated data and generic gene labels rather than any real findings
from that manuscript.

**This script uses simulated, randomly generated expression data, not
real patient data, real gene identities, or any third-party code.** It
was written from scratch to illustrate the method only.

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

This demo illustrates the type of pan-cancer meta-analysis screening
approach used in biomarker discovery research, including a
collaborative manuscript I co-authored that is currently under review
at a peer-reviewed journal. Specific findings from that study are not
included here, as it has not yet been published.
