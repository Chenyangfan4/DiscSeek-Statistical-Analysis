# Study dataset construction and statistical analysis conventions

## External data organization

The external CSV contains one selected examination per subject, represented by up to six cervical disc-level rows. The scripts use `ImageID` to identify an examination.

For **descriptive** analysis (`02_descriptive_analysis.R`), the H1/S1 Manual subset is represented both in its own group (n = 225) and as part of the complete Multi-center Manual group (n = 889). This overlapping representation is for group-specific summary statistics and concordance assessments, not independent comparisons.

For **scan angle** and **classifier** analysis, records are not duplicated. `ScanType == "Manual"` identifies the complete manual cohort; the separate `ScanGroup` names indicate H1/S1 Manual and the non-H1/S1 manual subset. The group code `Manual` in modeling refers only to the non-H1/S1 manual subset, not all 889 examinations.

## Diagnosis-level statistics

Kappa, agreement rate, Bowker tests, and kappa comparisons are computed at the examination level (one row per examination, using the C2/C3 record or equivalently unique group/examination pairs). Agreement CI uses the uncorrected score interval from `stats::prop.test(correct = FALSE)`.

## Scan angle offset

The confirmed multicenter-manual mixed-effects formula is `Scan.angle.offset ~ ScanSite * IVDType + ScanMachine * IVDType + (1 | ImageID)` and its fixed effects are assessed by `car::Anova(type = "III")`. Cross-platform comparisons rely on mixed-effects model estimates; independent, platform-matched H1/S1 planning groups use `compare_means()` as implemented in the analysis script.

## Classification pipeline

The original experiment performs missForest imputation before patient-level train/test partitioning. This procedure is preserved exactly as a historical analysis choice; it may incorporate information from held-out cases into feature imputation and should not be represented as leakage-free preprocessing. SMOTE is applied to training data only; held-out test observations are not synthetically augmented. The shared examination-level partitions are saved as a local RDS for the R notebook. Bootstrap is configured for 1,000 resamples.

The staged XGBoost tuning process uses the output of preceding stages when selecting subsequent hyperparameters. Run the complete stages as documented in source comments before final evaluation. The notebook and XGBoost scripts use separate R package environments.
