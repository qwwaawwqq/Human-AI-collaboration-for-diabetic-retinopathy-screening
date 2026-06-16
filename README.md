# Human–AI collaboration for diabetic retinopathy screening

Bayesian bivariate diagnostic test-accuracy network meta-analysis (DTA-NMA) comparing
**autonomous AI**, **human readers**, and **human–AI collaboration** for diabetic
retinopathy screening from colour fundus photographs (40 studies).

Protocol: INPLASY202630009. Manuscript under review at *Telehealth and Digital Health*.

## Repository structure
```
.
├── data/      De-identified study-level 2×2 dataset (sheets: patient_rr, eye_rr) — no patient-level data
├── R/         Analysis & figure scripts (run from the repository root)
├── py/        Network-geometry figure (Python / matplotlib)
└── results/   Pre-computed fitted models (*.rds) and result tables (*.csv),
               provided so figures and headline numbers reproduce without re-fitting
```

| Script | Purpose |
|---|---|
| `R/DR_Pathway_NMA_RUN50k.R` | **Full pipeline** (4 chains × 50,000 iter / 10,000 burn-in): fits the models, generates Figures 2–4 (forest, HSROC, multivariate MR) and the supplementary convergence/funnel/MR plots, writes result CSVs and `results/*.rds` |
| `R/extract_stats.R` | Canonical pooled Se/Sp/Youden/SUCRA/P(best) and per-10,000 numbers, read from the fitted models |
| `R/nma_analysis.R` | Bivariate DTA-NMA: pooled Se/Sp, Youden, SUCRA/ranking, pairwise contrasts, unit-of-analysis sensitivity |
| `R/unit_of_analysis_v7.R` | Unit-of-analysis sensitivity across table-per-study specifications (Supplementary Table S14) |
| `R/leave_one_out_Li2024b.R` | Leave-one-out sensitivity omitting the dominant study (Li et al. 2024b, DeepDR-LLM) |
| `R/figure4_multivariate_MR.R` | Figure 4: combined patient + eye multivariate meta-regression forest plot |
| `R/metareg.R` | Univariate categorical meta-regression (26 covariates, patient & eye) |
| `R/prevalence.R` | Expected outcomes per 10,000 screened |
| `R/nma_with_lee.R` | Sensitivity analysis documenting the addition of Lee et al. 2021 |
| `py/make_network_fig.py` | Network-geometry figure (Supplementary Figure S15) |

## Reproduce
```bash
# from the repository root
Rscript R/DR_Pathway_NMA_RUN50k.R                          # full pipeline: models, Figures 2-4, CSVs
Rscript R/extract_stats.R                                  # headline Se/Sp/Youden/SUCRA (from results/*.rds)
Rscript R/nma_analysis.R patient_rr patient dataset_id 0   # patient-level NMA
Rscript R/nma_analysis.R eye_rr     eye     dataset_id 0   # eye-level NMA (primary)
Rscript R/unit_of_analysis_v7.R                            # Table S14 unit-of-analysis
Rscript R/leave_one_out_Li2024b.R                          # leave-one-out (drop Li 2024b)
Rscript R/metareg.R                                        # meta-regression
python  py/make_network_fig.py                             # network geometry (Figure S15)
```
The fitted models ship in `results/`, so `extract_stats.R` and the figure scripts run in
seconds without re-running the MCMC (which takes ~hours). Re-running
`DR_Pathway_NMA_RUN50k.R` will overwrite them.

## Requirements
- R (≥ 4.3) with `rjags`, `R2jags`, `coda`, `readxl`, `parallel`; and JAGS ≥ 4.3
  (https://mcmc-jags.sourceforge.io). `rjags` locates JAGS automatically.
- Python (≥ 3.9) with `matplotlib` (for `py/make_network_fig.py` only).

## Data
Aggregate, study-level 2×2 counts and study characteristics extracted from published
reports only — no individual-patient data.

## Citation & license
See `CITATION.cff`. Released under the MIT License (`LICENSE`).
