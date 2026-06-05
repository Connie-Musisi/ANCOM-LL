# ANCOM-LL: A Log-Linear Variant of ANCOM-BC

ANCOM-LL is a log-link alternative to [ANCOM-BC](https://doi.org/10.1038/s41467-020-17041-7) for differential abundance (DA) testing in microbiome studies. It retains ANCOM-BC's sampling-fraction bias correction but replaces the log-count linear model with a direct model for `log E(N)`, avoiding pseudo-count dependence and targeting a more interpretable mean-scale log-fold change.

This work is part of a PhD thesis on adaptive methods for differential abundance analysis at Hasselt University (UHasselt).

---

## Motivation

ANCOM-BC estimates group differences from mean log-counts, targeting `E[log N]`. Under overdispersion and sparsity, `E[log N] ≠ log E[N]`, so estimated effects can be biased relative to the mean-scale LFC. ANCOM-LL models the count mean directly, which:

- avoids adding a pseudo-count before log-transformation,
- targets `log E(N)` more directly,
- allows variance to be estimated using procedures robust to mean-variance misspecification.

---

## Model

For taxon *j* in sample *i*, with group indicator *X_i* and unknown sampling fraction *c_i*:

```
E(N_ij) = c_i · θ_j · exp(X_i · α_j)
```

Taking logs:

```
log E(N_ij) = S_i + τ_j + X_i · α_j
```

where `S_i = log(c_i)` is the sampling-fraction offset and `α_j` is the LFC of interest. Parameters are estimated iteratively via a closed-form M-estimating procedure (not a full GLM) that alternates between estimating sampling fractions `D_k` and taxon baselines `μ_j`.

Because the model is not identified without a constraint, a data-driven **reference set** *J* is selected from taxa whose pairwise contrast variances fall in the lowest 10th percentile. The bias-corrected contrast `α_j − γ̂` is then used for inference, where `γ̂ = mean(μ1_J) − mean(μ0_J)`.

---

## Implementation

The pipeline has three components:

### 1. Parameter estimation — `est.par2(db)`

Iteratively estimates taxon baselines (`mu0`, `mu1`) and sampling fractions (`D0`, `D1`) from a long-format data frame with columns `O` (count), `group`, `subject`, and `taxon`.

### 2. Reference set selection — `ref(db)`

Calls `est.par2`, computes per-taxon pairwise contrast variances `V_j`, selects the bottom 10th percentile as the reference set *J*, and returns the bias correction term `gamma_hatdiff`.

### 3. Variance estimation — two options

**Ordinary bootstrap (`ANCOMLL`):** Resamples subjects with replacement within each group across `B` bootstrap replicates. The full `ref()` pipeline is re-run on each bootstrap dataset and the empirical variance of bias-corrected contrasts is used as the standard error.

**Wild bootstrap + smoothed variance (`WildBootstrap_SmoothedVariance`):** An earlier formulation operating on `phyloseq` objects. The wild bootstrap uses a polynomial weighting function optimised to match the distribution of estimating-equation residuals. The smoothed variance estimator pools the mean-variance relationship across taxa via a Poisson GLM, then uses Monte Carlo simulation from a bivariate normal to propagate uncertainty through the log transformation.

---

## Simulation results

Evaluated on 100 negative-binomial datasets (250 taxa, 100 samples per group, LFC = 1, 10% DA taxa, ~28% zeros):

| Method        | FDR   | Sensitivity |
|---------------|-------|-------------|
| ANCOM-BC      | 0.071 | 0.884       |
| ANCOM-LL SV   | 0.054 | 0.212       |
| ANCOM-LL WB   | 0.127 | 0.653       |

The nominal FDR level is 0.05. ANCOM-LL SV achieves FDR control but at the cost of very low sensitivity. ANCOM-LL WB recovers sensitivity but exceeds the nominal FDR. Both variance estimators systematically underestimate the true variance of the bias-corrected contrast — the wild bootstrap more severely — confirming that variance estimation rather than the link function is the dominant bottleneck.

---

## Repository structure

```
.
├── R/
│   ├── ANCOMLL-functions.R       # Core pipeline: sim.data, est.par2, ref, ANCOMLL
│   └── WildBootstrap_SmoothedVariance.R      # Phyloseq-based pipeline with wild bootstrap + smoothed variance
├── figures/
│   ├── sensitivity_fdr.R         # Plotting script for FDR/sensitivity bar charts
│   ├── FDR_comparison.png
│   └── sensitivity_fdr.png
└── README.md
```

---

## Dependencies

```r
install.packages(c("parallel", "doSNOW", "foreach", "MASS", "mvtnorm",
                   "ggplot2", "patchwork", "tidyr", "dplyr", "phyloseq"))
```

`mvtnorm` is required for `rmvnorm()` used in the smoothed-variance estimator.

---

## Quick start

```r
source("R/ANCOMLL-functions.R")

# Simulate data
dat <- sim.data(n = 100, n.taxa = 250, FC = c(1, 3, 2),
                phi = 0.5, SF0 = 0.8, SF1 = 1.2)

# Run ANCOM-LL with ordinary bootstrap (B replicates)
res <- ANCOMLL(dat$db, B = 200)

# DA taxa at 5% FDR
which(res$adjusted_p < 0.05)
```

For the phyloseq-based pipeline (wild bootstrap or smoothed variance):

```r
source("R/Version_09_03_2023.R")

# physeq must have a "DE.ind" column in tax_table and a group variable in sample_data
res_wb  <- run.scenario(physeq, V.method = "none", var.method = "wild",     B = 100)
res_sv  <- run.scenario(physeq, V.method = "none", var.method = "smoothed", B = 100)
```

---

## Citation

> Lin, H., Peddada, S.D. Analysis of compositions of microbiomes with bias correction. Nat Commun 11, 3514 (2020). https://doi.org/10.1038/s41467-020-17041-7
