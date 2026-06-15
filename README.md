# Human–AI collaboration for diabetic retinopathy screening

Bayesian bivariate diagnostic test-accuracy network meta-analysis (DTA-NMA) comparing
**autonomous AI**, **human readers**, and **human–AI collaboration** for diabetic
retinopathy screening from colour fundus photographs (40 studies).

Protocol: INPLASY202630009. Manuscript under review at *Telehealth and Digital Health*.

## Repository structure
```
.
├── data/      De-identified study-level 2×2 dataset (sheets: patient_rr, eye_rr) — no patient-level data
├── R/         Analysis scripts (run from the repository root)
└── results/   Output directory (created on run)
```

| Script | Purpose |
|---|---|
| `R/nma_analysis.R` | Bivariate DTA-NMA: pooled Se/Sp, Youden, SUCRA/ranking, pairwise contrasts, unit-of-analysis sensitivity |
| `R/metareg.R` | Univariate categorical meta-regression (26 covariates, patient & eye) |
| `R/prevalence.R` | Expected outcomes per 10,000 screened |
| `R/nma_with_lee.R` | Sensitivity analysis documenting the addition of Lee et al. 2021 |

## Reproduce
```bash
# from the repository root
Rscript R/nma_analysis.R patient_rr patient dataset_id 0   # patient-level
Rscript R/nma_analysis.R eye_rr     eye     dataset_id 0   # eye-level (primary)
Rscript R/metareg.R                                        # meta-regression
```
Outputs are written to `results/`.

## Requirements
R (≥ 4.3) with `rjags`, `coda`, `readxl`, `parallel`; and JAGS ≥ 4.3
(https://mcmc-jags.sourceforge.io). `rjags` locates JAGS automatically.

## Data
Aggregate, study-level 2×2 counts and study characteristics extracted from published
reports only — no individual-patient data.

## Citation & license
See `CITATION.cff`. Released under the MIT License (`LICENSE`).
