# DiscSeek | External study data cleaning
# Input: data/external/total_ivd.csv. Output: results/processed/processed_data.rds.

source("R/00_config.R")
library(tidyverse)
file_dir <- file.path(INPUT_DIR, "external", "total_ivd.csv")

# 1. Read external disc-level data ----
# Load the external prospective validation data.
df <- read.csv(file_dir) 

df <- df%>%
  mutate(across(c(ScanType, Clinical.diagnosis, IVDType, Sex, ScanMachine, ScanSite, Radiology.diagnosis), as.factor)) %>%
  mutate(across(OUTLIER_CONFIG$discrete_params, as.factor))

# 2. Screen scan angle offsets ----
# Screen records by the configured scan angle offset criterion.
step1_data <- filter_data(df)

# 3. Flag anatomical measurement outliers ----
# Flag outlying anatomical measurements.
step2_data <- step1_data %>% mark_all_outliers()

# 4. Apply outlier masking ----
# Replace outlying anatomical measurements with missing values.
clean_data <- step2_data %>% remove_all_outliers()


# 5. Export the cleaned study dataset ----
# Save the cleaned disc-level analysis data.
readr::write_rds(clean_data %>% dplyr::select(-dplyr::any_of(c("PatientName", "PatientID"))),
                 file.path(CLEAN_DATA_DIR, "processed_data.rds"))
readr::write_csv(clean_data %>% dplyr::select(-dplyr::any_of(c("PatientName", "PatientID"))),
                 file.path(CLEAN_DATA_DIR, "processed_data.csv"))
