# ANCOM-LL: A Log-Linear Variant of ANCOM-BC for Differential Abundance Analysis

## Overview

ANCOM-LL is a log-link alternative to [ANCOM-BC](https://doi.org/10.1038/s41467-020-17041-7) for differential abundance (DA) testing in microbiome studies. It retains ANCOM-BC's bias-correction rationale — estimating and removing a sample-specific sampling-fraction term before testing — but replaces the log-count linear model with a quasi-Poisson GLM with a log link. This means the model targets `log E(N)` rather than `E(log N)`, which is a more interpretable estimand under overdispersion and sparsity, and avoids the need for a pseudo-count when zeros are present.

The method was developed as part of a PhD thesis on adaptive methods for differential abundance analysis at Hasselt University (UHasselt), under the supervision of Olivier Thas.

## Motivation

ANCOM-BC estimates log-fold changes (LFCs) from the difference of mean log-counts across groups. When data are sparse or overdispersed, `E[log N] ≠ log E[N]`, so the estimated effect size can be biased relative to the mean-scale LFC. ANCOM-LL addresses this by modelling the count mean directly via a log-link GLM. Additionally, it avoids having to add a pseudo-count to handle zero observations, which can introduce arbitrary bias for rare taxa.

## Model

For taxon *j* in sample *i*, let *N_ij* be the observed count, *c_i* the unknown sampling fraction, *X_i* the binary group indicator, and *α_j* the log-fold change of interest. The mean model is:

```
log E(N_ij) = S_i + τ_j + X_i * α_j
```

where `S_i = log(c_i)` is the sampling-fraction offset and `τ_j` is the taxon baseline. Parameters are estimated via M-estimating equations under a quasi-Poisson working model, which is consistent under mean-model misspecification.

Because this model is not identified without a reference constraint, a data-driven reference set *J* is selected — taxa in the lowest 10th percentile of the distribution of pairwise contrast variances. The bias correction term *γ̂* is then estimated as the average estimated abundance difference between groups within *J*, and the bias-corrected contrast `α_j - γ̂` is used for inference.

## Variance Estimation

Variance estimation is the main methodological challenge. Three strategies were evaluated:

- **Sandwich (Huber-White)**: per-taxon robust variance, consistent under misspecification but unstable at low depth.
- **Smoothed mean-variance sandwich (SV)**: pools the mean-variance relationship across taxa via LOESS before substituting into the sandwich formula, borrowing strength for low-abundance taxa.
- **Wild bootstrap (WB)**: perturbs counts with Rademacher or Mammen multipliers and re-estimates the full pipeline (including reference-set selection) in each replicate, propagating uncertainty from the bias-correction step.

## Results Summary

Evaluated on negative-binomial simulations (250 taxa, 100 samples, LFC = 1, 10% DA taxa, ~28% zeros):

- **ANCOM-BC**: FDR near nominal (5%), competitive sensitivity — used as the reference benchmark.
- **ANCOM-LL SV**: FDR near nominal but substantially reduced sensitivity (conservative behaviour).
- **ANCOM-LL WB**: Recovers some sensitivity but with inflated FDR (loss of calibration).

Both variance estimators systematically underestimate the true variance of the bias-corrected contrast, the wild bootstrap more severely. This is the primary bottleneck: the variance estimation problem on the count scale is harder than on the log-count scale, and neither estimator adequately propagates the uncertainty arising from reference-set selection.

The conclusion is that switching from a log-count to a log-link model does not produce material gains once the bias correction and reference frame are accounted for. The reference frame — not the link function — is the dominant driver of performance in both ANCOM-BC and ANCOM-LL.

## Repository Structure

```
.
├── R/
│   ├── ancomll_fit.R          # Main fitting function (quasi-Poisson GLM + bias correction)
│   ├── ancomll_variance.R     # Sandwich, smoothed-variance, and wild-bootstrap estimators
│   └── reference_selection.R  # Data-driven reference set selection
├── simulations/
│   ├── nb_simulation.R        # Negative-binomial simulation study
│   └── eval_ancomll.R         # Evaluation metrics (FDR, sensitivity)
├── vignettes/
│   └── ancomll_demo.Rmd       # Reproducible demo on the Dietswap dataset
└── README.md
```

## Dependencies

```r
install.packages(c("phyloseq", "MASS", "sandwich", "lmtest", "microbiome"))
```

The `microbiome` package (Bioconductor) is used for the Dietswap example dataset.

## Quick Start

```r
source("R/ancomll_fit.R")

# count_mat: taxa × samples matrix of non-negative integers
# group:     binary group vector (0/1), length = ncol(count_mat)

result <- ancomll(count_mat, group, variance = "sv", alpha = 0.05)
head(result$results)  # taxon-level LFCs, standard errors, q-values
```

## Citation

This work is part of:

> Musisi, C., Thas, O., Jaspers, S., Kodalci, L. and Babiera, J. (2026). An Adaptive Test for Differential Abundance in Microbiome Studies. *Submitted to PLOS Computational Biology* (PCOMPBIOL-D-25-01876).

## Author

Connie Musisi  
PhD candidate, Hasselt University (UHasselt)  
GitHub: [@Connie-Musisi](https://github.com/Connie-Musisi)
