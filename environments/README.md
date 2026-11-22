# Analysis environments

The R scripts were run with Windows 10 and R 4.4.2; the sensitivity-analysis notebook used Ubuntu 24.04.1, Conda, and R 4.4.2. The supplied files are historical environment records, not tested lockfiles.

The Conda package manifest reports `r-caret=6.0_94`, whereas the R session loaded `caret 7.0-1` from the same Conda prefix (confirmed using `packageVersion()` and `find.package()`). The runtime package version should be treated as authoritative for the original notebook. Do not assume recreating this YAML alone will reproduce the runtime exactly. The Windows R script environment used `xgboost 1.7.9.1`; the Jupyter environment loaded `xgboost 3.2.0.1`.

Local home-directory paths have been redacted from the session logs and replaced in the YAML. The environment manifests may include unrelated software because they describe complete historical environments.
