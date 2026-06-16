# Human–AI collaboration for diabetic retinopathy screening

Bayesian bivariate diagnostic test-accuracy network meta-analysis (DTA-NMA) comparing
**autonomous AI**, **human readers**, and **human–AI collaboration** for diabetic
retinopathy screening from colour fundus photographs (40 studies).

Protocol: INPLASY202630009. Manuscript under review at *Telehealth and Digital Health*.

## Repository structure
```
.
├── data/      De-identified study-level 2×2 dataset (sheets: patient_rr, eye_rr) — no patient-level data
├── R/         dr_nma_analysis.R — single, end-to-end analysis script (run from the repository root)
├── py/        make_network_fig.py — network-geometry figure (Supplementary Figure S15)
└── results/   Fitted models (*.rds) from the full run + sensitivity logs; all other
               outputs (CSVs, Figures 2–4) are written here when the script runs
```

## One script does everything: `R/dr_nma_analysis.R`
Run from the repository root:
```bash
Rscript R/dr_nma_analysis.R                  # full analysis (4 chains × 50,000 iter / 10,000 burn-in)
DR_QUICK=1 Rscript R/dr_nma_analysis.R       # fast smoke-test (tiny MCMC, ~1 min) — for checking the pipeline
python  py/make_network_fig.py               # network-geometry figure (Figure S15)
```
It reproduces, writing everything to `results/`:

| Section | Outputs |
|---|---|
| Primary NMA (patient & eye) | pooled Se/Sp/Youden, SUCRA + P(best) ranking, pairwise contrasts, heterogeneity, convergence → `*_per_arm.csv`, `*_pairwise.csv`, `*_ranking.csv`, `*_het.csv`, `*_conv.csv`, `*_draws.rds` |
| Sensitivity (eye) | unit-of-analysis (base / by-publication / one-per-dataset×arm / one-per-publication×arm) **and** leave-one-out dropping the dominant study (Li 2024b) → `eye_sensitivity.csv` |
| Expected outcomes | per 10,000 screened across a prevalence grid → `*_prevalence.csv` |
| Meta-regression | univariate (26 covariates) and multivariate → `metareg_*.csv`, `metareg_*_multivar.csv` |
| Figures | Figure 2 (forest), Figure 3 (HSROC), Figure 4 (multivariate meta-regression) → `Figure2/3/4_*.tiff` |

Selection rules match the manuscript: the unit-of-analysis sensitivity keeps the **largest**
2×2 table per group×arm; collaboration remains highest-ranked across all specifications and
under leave-one-out (SUCRA ≈ 0.89–0.94), though its advantage over AI alone is not statistically
significant. `results/` ships with the fitted posterior models (`AI_NMA_model.rds`,
`eye_rr_NMA_model.rds`) and the unit-of-analysis / leave-one-out result logs from the full run.

## Requirements
- R (≥ 4.3) with `rjags`, `coda`, `readxl`, `parallel`; and JAGS ≥ 4.3
  (https://mcmc-jags.sourceforge.io). `rjags` locates JAGS automatically.
- Python (≥ 3.9) with `matplotlib` (for `py/make_network_fig.py` only).

## Data
Aggregate, study-level 2×2 counts and study characteristics extracted from published
reports only — no individual-patient data.

## Citation & license
See `CITATION.cff`. Released under the MIT License (`LICENSE`).
