# FDA Oncology Drug Approval Analysis

Statistical Analysis of Priority Review Effects on FDA Drug Approval Timelines 

## Honors Project
~ *Experimental Design for Machine Learning Statistics (STA 4157) Project*

By: Basir Abdul Samad

## Overview

This repository contains a comprehensive statistical analysis investigating whether oncology drugs receive differential reductions in FDA approval timelines from Priority Review designation compared to non-oncology therapeutic areas. The analysis reveals Simpson's Paradox in the therapeutic area main effect when controlling for regulatory era.

## Key Findings

- **Simpson's Paradox Detected**: The therapeutic area main effect decreased from F=21.22 (unblocked model) to F=18.70 (era-blocked model), demonstrating that part of the apparent oncology speed advantage was confounded with temporal trends
- **No Interaction Effect**: Analysis found no statistically significant interaction between therapeutic area and review designation across all models (Model 1: F=0.08, p=0.776; Model 2: F=0.32, p=0.572; Model 3: F=2.43, p=0.119), confirming Priority Review provides consistent temporal benefit across therapeutic areas
- **Mechanism Identified**: The perceived speed advantage for oncology drugs reflects strategic utilization of the Accelerated Approval pathway (4.9× higher usage: 30.4% vs 6.2%) rather than systemic FDA bias
- **Strong Main Effects**: Review designation (F=77.52, $p < 0.001$), Regulatory era (F=88.64, $p < 0.001$), Therapeutic area (F=18.70, $p < 0.001$)

## Repository Structure

```
fda-oncology-approval-analysis/
├── docs/                          # Final documentation
│   ├── FDA_Analysis_Guide_FINAL.tex
│   ├── FDA_Analysis_Guide_FINAL.pdf
│   └── fda_references.bib
├── data/                          # FDA dataset
│   └── final_for_posting_compilation_of_cder_nme_and_new_biologic_approvals_1985-2024.csv
├── visualizations/                # Publication-quality figures (22 figures)
│   ├── fig1_interaction_plot.png
│   ├── fig2_era_trends_plot.png
│   ├── fig3_permutation_histograms.png
│   ├── fig4_cell_means_barchart.png
│   ├── fig5_effect_sizes_barchart.png
│   └── ... (14 more figures)
├── results/                       # Analysis output files
│   ├── fda_analysis_clean.csv
│   ├── fda_analysis_clean_classified.csv
│   ├── model_comparison.csv
│   ├── cell_means.csv
│   ├── permutation_test_results.csv
│   └── ... (10+ CSV files)
├── scripts/                       # R analysis pipeline
│   ├── config.R
│   ├── utils.R
│   ├── phase1_data_loading.R
│   ├── phase2_classification.R
│   ├── phase3_missing_data.R
│   ├── phase4_balance_validation.R
│   ├── phase4b_expedited_programs.R
│   ├── phase5_assumption_testing.R
│   ├── phase6_anova_modeling.R
│   ├── phase6b_ancova_modeling.R
│   ├── phase7_era_stratification.R
│   ├── phase8_sensitivity_analysis.R
│   └── phase8b_sensitivity_uncertain.R
└── README.md                      # This file
```

## Analysis Pipeline

The analysis follows an 11-phase R pipeline:

1. **Phase 1**: Data loading and cleaning (n=1,335 drugs)
2. **Phase 2**: Therapeutic area classification (53 oncology keywords)
3. **Phase 3**: Missing data analysis (MCAR assessment)
4. **Phase 4**: Balance and structural validation (n=1,038 analysis-ready)
5. **Phase 4b**: Expedited programs analysis (Accelerated Approval, Fast Track, Breakthrough, Orphan)
6. **Phase 5**: Assumption testing (log transformation validation)
7. **Phase 6**: ANOVA modeling (Type III SS, era blocking)
8. **Phase 6b**: ANCOVA modeling (expedited programs as covariates)
9. **Phase 7**: Era stratification (confirming no interaction across all eras)
10. **Phase 8**: Sensitivity analysis (high-confidence subset)
11. **Phase 8b**: Sensitivity analysis on uncertain classifications

## Statistical Methods

- **Two-way factorial ANOVA** with Type III Sums of Squares
- **Era blocking** to control for temporal trends (Pre-PDUFA, Early-PDUFA, Mid-PDUFA, Post-FDASIA)
- **Log transformation** for variance stabilization
- **Permutation testing** (10,000 iterations) for assumption-free inference
- **ANCOVA** with expedited program covariates
- **Sensitivity analysis** on high-confidence and uncertain classifications

## Requirements

### R Packages
- tidyverse
- car (Type III SS)
- emmeans
- ggplot2

### LaTeX
- pdflatex
- BibTeX
- Packages: natbib, amsmath, booktabs, graphicx, hyperref

## Usage

### Running the Analysis

```r
# Run all phases sequentially
source("scripts/phase1_data_loading.R")
source("scripts/phase2_classification.R")
source("scripts/phase3_missing_data.R")
source("scripts/phase4_balance_validation.R")
source("scripts/phase4b_expedited_programs.R")
source("scripts/phase5_assumption_testing.R")
source("scripts/phase6_anova_modeling.R")
source("scripts/phase6b_ancova_modeling.R")
source("scripts/phase7_era_stratification.R")
source("scripts/phase8_sensitivity_analysis.R")
source("scripts/phase8b_sensitivity_uncertain.R")
```

### Compiling the Document

```bash
cd docs/
pdflatex FDA_Analysis_Guide_FINAL.tex
bibtex FDA_Analysis_Guide_FINAL
pdflatex FDA_Analysis_Guide_FINAL.tex
pdflatex FDA_Analysis_Guide_FINAL.tex
```

## Data Source

U.S. Food and Drug Administration (2024). *Compilation of CDER new molecular entity (NME) drug and new biologic approvals*.
https://www.fda.gov/drugs/drug-approvals-and-databases/compilation-cder-new-molecular-entity-nme-drug-and-new-biologic-approvals

## Results Summary

| Analysis | Sample Size | Key Finding |
|----------|-------------|-------------|
| Pooled (Raw, no blocking) | n=1,038 | No interaction (F=0.08, p=0.776), Strong main effects |
| Pooled (Raw, era-blocked) | n=1,038 | No interaction (F=0.32, p=0.572), Simpson's Paradox in main effect |
| Pooled (Log, era-blocked) | n=1,038 | No interaction (F=2.43, p=0.119), Preferred model |
| Era-stratified | Per era | No interaction in any era (all p>0.10) |
| ANCOVA (Expedited Programs) | n=1,038 | Accelerated Approval main driver (30.4% vs 6.2%) |
| Sensitivity (High-confidence) | Subset | Findings robust to classification uncertainty |

## Key Visualizations

- **Figure 1**: Interaction plot showing differential Priority Review effects
- **Figure 2**: Temporal trends across regulatory eras
- **Figure 3**: Permutation test distribution (10,000 iterations)
- **Figure 4**: Cell means with 95% confidence intervals
- **Figure 5**: Effect sizes comparison (η²ₚ)
- **Figure 6**: Residual diagnostics for model validation
- **Additional figures**: Era-stratified analyses, ANCOVA results, sensitivity analyses

## Citation

```bibtex
@unpublished{samad2025fda,
  author = {Abdul Samad, Basir},
  title = {Statistical Analysis of Priority Review Effects on FDA Drug Approval Timelines: Evidence of Simpson's Paradox in Oncology Drug Regulation},
  year = {2025},
  note = {Honors Research Project}
}
```

## License

This project is for academic purposes. Data source: U.S. Food and Drug Administration (public domain).

---

**Analysis Date**: December 10, 2025
