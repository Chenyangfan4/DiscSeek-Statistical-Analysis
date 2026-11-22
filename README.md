# DiscSeek: statistical analysis and plotting code

This repository contains the non-proprietary R analysis scripts and an R-kernel Jupyter notebook supporting the DiscSeek cervical spine MRI scan-planning study. It covers independent internal technical validation and external prospective validation with multicenter, multivendor **manual** planning benchmarks. DiscSeek was deployed only on the platform-matched H1/S1 system. MRI scanner integration, V-Net training source, and model weights are not included.

## Inputs and privacy

Controlled-access MRI/clinical data are not distributed. Prepare the three input CSV files with exact column names described in [`docs/data_dictionary.md`](docs/data_dictionary.md). Each external CSV row is a **disc-level observation**, not a separate subject; examinations contribute multiple disc-level rows. Patient identifiers and derived row-level data must not be deposited publicly.

- `data/external/total_ivd.csv`: external prospective validation, disc-level records.
- `data/internal/internal_planning_error.csv`: internal retrospective technical validation, one test examination per row.
- `data/training/loss_metric_curves.csv`: training/validation logs from the parallel V-Nets.

## Environments

Main R scripts: R 4.4.2 under Windows 10/RStudio. Sensitivity-analysis notebook: R 4.4.2, Ubuntu 24.04.1, Conda/IRkernel. Environment snapshots and the documented caret manifest/runtime difference are in [`environments/`](environments/README.md). These records have not been independently rebuilt.

## Execution order

Open `DiscSeek-Analysis.Rproj` at the repository root. The scripts assume the R working directory is the repository root. All scripts have navigable RStudio `# Section ----` headings.

**Internal validation (separate branch):**

```r
source("R/internal/internal_analysis.R")       # Table 1, Figures S2–S3
source("R/internal/training_visualization.R") # Figure S4
```

**External validation:**

```r
source("R/01_preprocessing.R")             # clean data
source("R/02_descriptive_analysis.R")      # Tables 2, S1, S3–S5; Figure 2k–l
source("R/03_scan_angle_analysis.R")       # Table 2, S2; Figure 2i–j
source("R/04_modeling_data_preparation.R") # Table S9, shared model partitions
```

Then run `R/05a_xgboost_analysis.R` for Table 3, Figure 3 and Table S10, and run `notebooks/05b_classifier_sensitivity_analysis.ipynb` using the separate Conda R kernel for Table S6. Both depend on the partition output from `04_modeling_data_preparation.R`; the notebook reads `results/intermediate/resampling_data.rds`. The optional final cross-classifier visualization in 05a should be run only **after** the notebook exports its results. The staged tuning section requires retaining the chosen hyperparameters between stages.

Detailed study group definitions and limitations are in [`docs/data_processing.md`](docs/data_processing.md). An item-by-item mapping covering every Main and Supplementary Figure/Table is in [`docs/analysis_index.md`](docs/analysis_index.md).

## Output and reproduction scope

Generated tables, figures and intermediate analysis objects are written under `results/`. The scripts prepare the underlying statistical outputs and panels, not the final journal graphic assembly. Without the controlled-access input datasets, numeric findings cannot be independently reproduced from the public archive alone. The archive has been statically reviewed but not run end-to-end with source data.

## Citation and rights

Intended resource type for Zenodo: **Software**. Before publication, the rights holder must approve the license and creator metadata and add a version number and Zenodo DOI. A license file is intentionally not invented. Cite the deposited release DOI once available.
# HelloWorld
