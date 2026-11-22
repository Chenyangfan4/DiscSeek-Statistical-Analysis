# DiscSeek | External prospective validation: modeling data preparation
# Outputs: Supplementary Table S9; shared examination-level train/test data for 05a and 05b.
# This script retains the study's original pre-split missForest imputation protocol.

source("R/00_config.R")
library(tidyverse)
library(missForest)
library(caret)
library(DMwR)

# 1. Load and prepare external validation data ----
# Load preprocessed disc-level observations.
data <- read_rds(file.path(CLEAN_DATA_DIR, "processed_data.rds"))

# Convert categorical variables and set factor reference levels.
data <- data %>%
  mutate(across(c("ScanType", "ScanMachine", "ScanSite"), as.factor)) %>%
  mutate(across(contains("diagnosis"), as.factor)) %>%
  mutate(across(all_of(OUTLIER_CONFIG$sequence_params), ~ factor(.x, levels = sort(unique(.x)), ordered = TRUE))) %>%
  mutate(across(all_of(OUTLIER_CONFIG$discrete_params), as.factor))

data <- data %>%  
  mutate(
    IVDType = factor(data$IVDType, levels = sort(unique(data$IVDType))),
    Sex = factor(Sex)
  )

data$ScanMachine <- relevel(data$ScanMachine, ref="Philips")
data$ScanSite <- relevel(data$ScanSite, ref="QQHE")
data$ScanType<- relevel(data$ScanType, ref = "Manual")
data$Clinical.diagnosis <- relevel(factor(data$Clinical.diagnosis), ref = "ACS")
data$Radiology.diagnosis <- relevel(factor(data$Radiology.diagnosis), ref = "ACS")
data$IVDType <- relevel(data$IVDType, ref = "C2C3")

data_cleaned <- data %>%
  unite("ScanGroup", ScanSite, ScanMachine, ScanType, remove = FALSE) %>%
  mutate(
    ScanGroup = factor(ScanGroup),
    Age = (Age - min(data$Age)) / (max(data$Age) - min(data$Age)),
    Imaging.quality.score = factor(Imaging.quality.score, levels = sort(unique(data$Imaging.quality.score)))
  )

# 2. Construct patient-level model inputs and partition data ----
params <- c(OUTLIER_CONFIG$continuous_params, 
            OUTLIER_CONFIG$sequence_params, 
            OUTLIER_CONFIG$discrete_params)  # Anatomical predictors: continuous, ordered, and nominal measurement variables.
diag_var <- "Clinical.diagnosis"  # Clinical diagnosis is the classification target.
patient_var <- c("ImageID")  # Use ImageID to identify the selected examination.
scan_groups <- levels(data_cleaned$ScanGroup)  # Available acquisition-planning groups.
other_params <- c("Sex", "Age","IVDType")


# Impute missing anatomical measurements using missForest before the data split (original analysis workflow).
imputed_data <- missForest(data_cleaned[,params])$ximp
data_complete <- imputed_data
# Combine imputed anatomical measurements with examination-level labels and covariates.
data_complete <- data_complete %>% bind_cols(data_cleaned[,c(diag_var, patient_var, "ScanGroup")])
data_complete <- data_complete %>% bind_cols(data_cleaned[,other_params])
patient_data <- data_complete %>% pivot_wider(names_from = "IVDType", values_from = params)
# Configure the classifier family.
model_type <- "xgboost"
# Set up patient-level train/test partitions.
patient_scan_groups <- c("DT_UIH_AI", "DT_UIH_Interactive", "DT_UIH_Manual", "Manual")
train_val_list <- list(training=NULL, validation=NULL)
for (scan_group in scan_groups) {
  df <- split_preserving_diagnosis_caret(patient_data %>% filter(str_detect(ScanGroup, scan_group)), 
                                         "ImageID", diag_var, ratio = 0.8)
  train_val_list$training <- bind_rows(train_val_list$training, df$training)
  train_val_list$validation <- bind_rows(train_val_list$validation, df$validation)
}

sample_df <- train_val_list$training
# 3. Balance the H1/S1 planning groups with SMOTE ----
# Balance clinical diagnosis categories within each H1/S1 planning group.
ordered_cols <- sample_df %>%
  dplyr::select(where(~ "ordered" %in% class(.))) %>%
  names()

level_info <- lapply(ordered_cols, function(col) {
  list(
    levels = levels(sample_df[[col]])
  )
})
names(level_info) <- ordered_cols
balanced_data <- sample_df %>% 
  filter(str_detect(ScanGroup, "UIH")) %>%
  group_by(ScanGroup) %>%
  group_modify(~ {
    set.seed(123)
    df <- .x
    df_fixed <- df %>% dplyr::select(-ImageID) %>% as.data.frame()
    for (col in ordered_cols) {
      df_fixed[[col]] <- factor(df_fixed[[col]], levels = level_info[[col]]$levels, ordered = FALSE)
    }
    df_size <- table(df_fixed$Clinical.diagnosis)
    majority_cls <- names(which.max(df_size))
    minority_cls <- names(which(df_size < max(df_size)))
    majority_df <- df_fixed %>% filter(Clinical.diagnosis == majority_cls)
    smote_result <- NULL
    ratio = 1.0  # Majority-class sampling ratio.
    for (cls in minority_cls) {
      minority_df <- df_fixed %>% filter(Clinical.diagnosis == cls)
      # print(unique(minority_df$Clinical.diagnosis))
      new_df <- bind_rows(minority_df, majority_df)
      new_df$Clinical.diagnosis <- factor(new_df$Clinical.diagnosis, levels = c(majority_cls, cls))
      smote_data <- SMOTE(
        Clinical.diagnosis ~ .,
        data = new_df,
        K = 5,  # Number of nearest neighbors for SMOTE.
        perc.over = round((nrow(majority_df) * ratio - nrow(minority_df)) / nrow(minority_df) * 100),  # Generate synthetic minority-class observations to approach the target size.
        perc.under = 100  # Set the undersampling argument (only minority-class outputs are retained).
      )
      smote_data <- union(smote_data %>% filter(Clinical.diagnosis == cls), minority_df)
      print("+++++++++++++++++++++++++++++++++++++")
      print(paste("major", nrow(majority_df)))
      print(paste("minor", nrow(minority_df)))
      print(paste("smote", nrow(smote_data)))
      smote_result <- bind_rows(smote_result, smote_data)
    }
    # smote_size <- round(mean(table(smote_result$Clinical.diagnosis)))
    # majority_df <- majority_df %>% sample_n(smote_size)
    bind_rows(smote_result, majority_df)
  }) %>% ungroup()

balanced_data$ScanGroup <- factor(balanced_data$ScanGroup, levels = unique(balanced_data$ScanGroup))

# Balance sampling across H1/S1 planning groups within each diagnosis class.
upper_ratio <- 1.4 # Additional synthetic sampling ratio.
balanced_data <- balanced_data %>% 
  group_by(Clinical.diagnosis) %>%
  group_modify(~ {
    set.seed(123)
    df <- .x
    df_size <- table(df$ScanGroup)
    majority_cls <- names(which.min(df_size))
    unique(df$ScanGroup)
    smote_result <- SMOTE(ScanGroup ~ ., data = df %>% as.data.frame(), perc.over = 250, perc.under = 100, K=5)
    smote_result <- smote_result %>% 
      filter(ScanGroup==majority_cls) %>%
      anti_join(df, by = names(smote_result)) %>%
      sample_n(round(upper_ratio * df_size[[majority_cls]])) %>%
      union(df)
  })


for (col in ordered_cols) {
  balanced_data[[col]] <- factor(
    balanced_data[[col]], 
    levels = level_info[[col]]$levels, 
    ordered = TRUE
  )
}
# 4. Construct an independent non-H1/S1 manual reference subset ----
# Select training observations from other hospital–vendor settings.
manual_reference <- balanced_data %>%
  filter(ScanGroup == "DT_UIH_Manual") %>%
  count(Clinical.diagnosis) %>%
  rename(target_n = n)  # Use H1/S1 Manual diagnosis counts as the sampling reference.

remaining_data <- sample_df %>%
  filter(!ScanGroup %in% unique(balanced_data$ScanGroup))  # Exclude the H1/S1 planning groups already balanced above.

available_groups <- remaining_data %>% distinct(ScanGroup) %>% pull(ScanGroup)

manual_counts <- balanced_data %>% 
  filter(ScanGroup == "DT_UIH_Manual") %>% 
  count(Clinical.diagnosis)

balanced_remaining <- manual_counts%>% 
  {set.seed(123); .} %>% 
  distinct(Clinical.diagnosis, n) %>% 
  pmap_dfr(~ {
    cls <- ..1
    target_n <- ..2
    print(paste(cls, target_n))
    
    datasets <- unique(remaining_data %>% 
                         filter(Clinical.diagnosis == cls) %>% 
                         pull(ScanGroup))
    per_dataset <- floor(target_n / length(datasets))
    
    sampled <- remaining_data %>%
      filter(Clinical.diagnosis == cls) %>%
      mutate(ScanGroup = as.character(ScanGroup)) %>%
      group_by(ScanGroup) %>%
      sample_n(size = min(per_dataset, n()), replace = FALSE) %>%
      ungroup()
    
    extra <- target_n - nrow(sampled)
    while (extra > 0) {
      extra_dataset <- anti_join(remaining_data, sampled, by=names(sampled)) %>%
        filter(Clinical.diagnosis == cls) %>%
        pull(ScanGroup) %>%
        unique()
      extra_sampled <- anti_join(remaining_data, sampled, by=names(sampled)) %>%
        filter(Clinical.diagnosis == cls) %>%
        mutate(ScanGroup = as.character(ScanGroup)) %>%
        group_by(ScanGroup) %>%
        sample_n(size = max(floor(extra / length(extra_dataset)), 1), replace = FALSE) %>%
        ungroup() %>% 
        sample_n(size = min(extra, n()), replace = FALSE)
      print(paste("total:", per_dataset, "extra:", extra, "sampled:", nrow(sampled), "extra_sampled:", nrow(extra_sampled)))
      sampled <- bind_rows(sampled, extra_sampled)
      extra <- target_n - nrow(sampled)
    }
    print(paste("total:", per_dataset, "sampled:", nrow(sampled)))
    return(sampled)
  })
  

balanced_remaining$ScanGroup <- "Manual"

balanced_data <- bind_rows(balanced_data, balanced_remaining) %>% arrange(ScanGroup, Clinical.diagnosis)
balanced_data$ScanGroup <- factor(balanced_data$ScanGroup) 
# balanced_data$ImageID <- 1:nrow(balanced_data)


# Assign new ImageID values to synthetic training observations.
inner_data <- semi_join(sample_df, balanced_data, by=setdiff(names(sample_df), c("ScanGroup", "ImageID"))) %>% mutate(ScanGroup=as.character(ScanGroup))
outer_data <- anti_join(balanced_data, sample_df, by=setdiff(names(sample_df), c("ScanGroup", "ImageID"))) %>% mutate(ScanGroup=as.character(ScanGroup))
outer_data$ImageID <- (max(patient_data$ImageID)+1):(max(patient_data$ImageID)+nrow(outer_data))
balanced_data <- bind_rows(inner_data, outer_data) %>% mutate(ScanGroup = ifelse(str_detect(ScanGroup, "UIH"), ScanGroup, "Manual"))
balanced_data$ScanGroup <- factor(balanced_data$ScanGroup, levels = unique(balanced_data$ScanGroup))

train_val_list$balanced_training <- balanced_data  # Apply resampling only to the training partition; keep held-out observations unmodified.

# Sample the non-H1/S1 manual test subset independently.
remaining_data <- train_val_list$validation %>%
  filter(!str_detect(ScanGroup, "UIH"))  # Select manual observations outside the H1/S1 platform.

manual_counts <- train_val_list$validation %>% 
  filter(ScanGroup == "DT_UIH_Manual") %>% 
  count(Clinical.diagnosis)

balanced_remaining <- manual_counts%>% 
  {set.seed(123); .} %>% 
  distinct(Clinical.diagnosis, n) %>% 
  pmap_dfr(~ {
    cls <- ..1
    target_n <- ..2
    print(paste(cls, target_n))
    
    datasets <- unique(remaining_data %>% 
                         filter(Clinical.diagnosis == cls) %>% 
                         pull(ScanGroup))
    per_dataset <- floor(target_n / length(datasets))
    
    sampled <- remaining_data %>%
      filter(Clinical.diagnosis == cls) %>%
      mutate(ScanGroup = as.character(ScanGroup)) %>%
      group_by(ScanGroup) %>%
      sample_n(size = min(per_dataset, n()), replace = FALSE) %>%
      ungroup()
    
    extra <- target_n - nrow(sampled)
    while (extra > 0) {
      extra_dataset <- anti_join(remaining_data, sampled, by=names(sampled)) %>%
        filter(Clinical.diagnosis == cls) %>%
        pull(ScanGroup) %>%
        unique()
      extra_sampled <- anti_join(remaining_data, sampled, by=names(sampled)) %>%
        filter(Clinical.diagnosis == cls) %>%
        mutate(ScanGroup = as.character(ScanGroup)) %>%
        group_by(ScanGroup) %>%
        sample_n(size = max(floor(extra / length(extra_dataset)), 1), replace = FALSE) %>%
        ungroup() %>% 
        sample_n(size = min(extra, n()), replace = FALSE)
      print(paste("total:", per_dataset, "extra:", extra, "sampled:", nrow(sampled), "extra_sampled:", nrow(extra_sampled)))
      sampled <- bind_rows(sampled, extra_sampled)
      extra <- target_n - nrow(sampled)
    }
    print(paste("total:", per_dataset, "sampled:", nrow(sampled)))
    return(sampled)
  })


balanced_remaining$ScanGroup <- "Manual"

train_val_list$balanced_validation <- train_val_list$validation %>% 
  filter(str_detect(ScanGroup,"UIH")) %>% 
  bind_rows(balanced_remaining) %>% 
  arrange(ScanGroup, Clinical.diagnosis) %>%
  mutate(ScanGroup = factor(ScanGroup))

# Save the shared resampling object for the R notebook.
readr::write_rds(train_val_list, file.path(INTERMEDIATE_DIR, "resampling_data.rds"))

# 5. Supplementary Table S9 | Post-resampling class distribution ----
total_data = bind_rows(
  train_val_list$balanced_training %>% 
    mutate(dataset="training"), 
  train_val_list$balanced_validation %>% 
    mutate(dataset="validation"))
resampling_summary <- total_data %>% 
  group_by(ScanGroup, Clinical.diagnosis, dataset) %>% 
  summarise(n=n(), .groups='drop') %>% 
  group_by(ScanGroup) %>% 
  mutate(total=sum(n), 
         percentage=100*n/total, 
         result = sprintf("%d(%.0f)%%", n, percentage)) %>% 
  ungroup() %>% 
  dplyr::select(-c(n, total, percentage)) %>%
  pivot_wider(names_from = ScanGroup, values_from = result) %>%
  arrange(dataset)


readr::write_csv(resampling_summary, file.path(MODELING_TABLE_DIR, "table_s9_resampling_distribution.csv"))

# 6. Export shared classifier inputs and feature metadata ----
readr::write_rds(
  list(train_val_list = train_val_list, patient_data = patient_data,
       patient_scan_groups = patient_scan_groups, diag_var = diag_var,
       patient_var = patient_var, model_type = model_type),
  file.path(INTERMEDIATE_DIR, "modeling_context.rds")
)
