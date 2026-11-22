# DiscSeek | Shared configuration and reusable analysis functions
# Use from repository root. Main scripts load their analysis-specific dependencies.

library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(ggpubr)
library(rstatix)
library(nortest)
library(lme4)
library(caret)
library(xgboost)
library(randomForest)
library(shapviz)
# Run the scripts from the repository root, with DiscSeek-Analysis.Rproj open.
PROJECT_DIR <- getwd()
INPUT_DIR <- file.path(PROJECT_DIR, "data")
RESULTS_DIR <- file.path(PROJECT_DIR, "results")
CLEAN_DATA_DIR <- file.path(RESULTS_DIR, "processed")
FIGURE_DIR <- file.path(RESULTS_DIR, "figures")
TABLE_DIR <- file.path(RESULTS_DIR, "tables")
INTERMEDIATE_DIR <- file.path(RESULTS_DIR, "intermediate")
MODELING_TABLE_DIR <- file.path(TABLE_DIR, "modeling_prediction")
COMPARISON_TABLE_DIR <- file.path(TABLE_DIR, "comparisons", "overall")
MODELING_FIGURE_DIR <- file.path(FIGURE_DIR, "modeling_prediction")
COMPARISON_FIGURE_DIR <- file.path(FIGURE_DIR, "comparisons", "overall")

# Create analysis output directories ----
for (output_dir in c(RESULTS_DIR, CLEAN_DATA_DIR, FIGURE_DIR, TABLE_DIR,
                     INTERMEDIATE_DIR, MODELING_TABLE_DIR, COMPARISON_TABLE_DIR,
                     MODELING_FIGURE_DIR, COMPARISON_FIGURE_DIR)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

# Outlier screening parameters
OUTLIER_CONFIG <- list(
  # Scan angle offset screening
  scan_offset = list(
    method = "iqr_image",
    max_image = 10,
    max_disc = 5,
    iqr_coef = 1.5
  ),
  
  # Anatomical parameter outlier screening
  param_outlier = list(
    methods = c("iqr", "zscore"),
    group_vars = c("IVDType", "Clinical.diagnosis"),
    z_thresh = 3,
    iqr_coef = 1.5
  ),
  
  # Ordered anatomical parameters
  sequence_params = c("Grade.of.IVD.protrusion",
                      "Grade.of.left.foraminal.stenosis",                
                      "Grade.of.right.foraminal.stenosis",
                      "Grade.of.spinal.stenosis"),
  
  # Nominal anatomical parameters
  discrete_params = c("Morphology.symmetry.of.bilateral.lateral.recesses"),
  
  # Continuous anatomical parameters
  continuous_params = c("Height.of.IVD..cm.",
                       "Height.of.left.lateral.recess..mm.",
                       "Height.of.right.lateral.recess..mm.",
                       "Sagittal.diameter.of.IVD..cm.",
                       "Sagittal.diameter.of.VB..cm.",
                       "Sagittal.diameter.of.spinal.canal..cm.",  
                       "Torg.ratio",
                       "Transverse.diameter.of.IVD..cm.",                  
                       "Transverse.diameter.of.spinal.canal..cm.")
  )


# Filter observations with excessive scan angle offsets using the configured screening method.
filter_data <- function(data) {
  cfg <- OUTLIER_CONFIG$scan_offset
  method <- cfg$method
  switch(method,
         "thresh_vert" = {
           data %>% filter(Scan.angle.offset <= cfg$max_disc)
         },
         "thresh_image" = {
           data %>% filter(Max.scan.angle.offset <= cfg$max_image)
         },
         "iqr_vert" = {
           data %>%
             group_by(IVDType) %>%
             mutate(upper = quantile(Scan.angle.offset, 0.75) + cfg$iqr_coef*IQR(Scan.angle.offset)) %>%
             filter(Scan.angle.offset <= upper) %>%
             ungroup()
         },
         "iqr_image" = {
           data %>%
             mutate(upper = quantile(Max.scan.angle.offset, 0.75) + cfg$iqr_coef*IQR(Max.scan.angle.offset)) %>%
             filter(Max.scan.angle.offset <= upper) %>%
             ungroup()
         }
  )
}

# Flag outlying anatomical measurements within disc-level and diagnosis strata. ----
mark_outliers <- function(data, param, group_vars = c("IVDType", "Clinical.diagnosis")) {
  # Test normality within each specified subgroup.
  norm_test <- data %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      p_value = shapiro.test(.data[[param]])$p.value,
      .groups = "drop"
    )
  
  # Optional: inspect subgroup normality results.
  # print(unique(norm_test))
  
  # Choose the outlier criterion based on subgroup normality and flag observations.
  data <- data %>%
    left_join(norm_test, by = group_vars) %>%
    group_by(across(all_of(group_vars))) %>%
    mutate(
      outlier_flag = if_else(
        p_value > 0.05,
        abs(scale(.data[[param]])) > 3,  # Z-score > 3
        .data[[param]] < quantile(.data[[param]], 0.25) - 1.5 * IQR(.data[[param]]) |
          .data[[param]] > quantile(.data[[param]], 0.75) + 1.5 * IQR(.data[[param]])
      )
    ) %>%
    ungroup() %>%
    select(-p_value) %>%
    rename_with(~ paste0(param, "_outlier"), "outlier_flag")
  
  # Optional: display the number of flagged values.
  # print(sum(data[[paste0(param, "_outlier")]], na.rm = TRUE))
  
  return(data)
}

# Apply outlier flags to all configured anatomical parameters. ----
mark_all_outliers <- function(data) {
  cfg <- OUTLIER_CONFIG$param_outlier
  for (param in OUTLIER_CONFIG$continuous_params) {
    data <- mark_outliers(data, param, cfg$group_vars)
  }
  return(data)
}

# Replace flagged anatomical values with missing values. ----
remove_all_outliers <- function(data) {
  # Identify all outlier indicator columns.
  outlier_cols <- names(data) %>% str_subset("_outlier$")
  
  # Set flagged values to NA for each anatomical parameter.
  for (col in outlier_cols) {
    param <- str_remove(col, "_outlier$")
    data <- data %>%
      mutate(
        !!param := if_else(.data[[col]], NA, .data[[param]])
      )
  }
  
  # Drop temporary outlier indicator columns.
  data %>%
    select(-all_of(outlier_cols))
}


# Select and run statistical comparisons for the requested groups. ----
perform_diff_test <- function(data, formula, group_var, is_norm, paired=FALSE, pairwise=TRUE) {
  # if (length(unique(data[[group_var]])) == 2) {
  # Optional: enable pairwise comparisons for two eligible groups.
  #   pairwise = TRUE
  # } 
  valid_group_vars <- data %>%
    group_by(!!sym(group_var)) %>%
    filter(n() >= 3) %>%
    pull(!!sym(group_var)) %>%
    unique()
  data <- data %>% filter(!!sym(group_var) %in% valid_group_vars)
  if (is_norm) {
    method <- ifelse(pairwise, "t.test", "anova")
    var.equal <- levene_test(data, formula = formula)$p > 0.05
  } else {
    method <- ifelse(pairwise, "wilcox.test", "kruskal.test")
    var.equal <- FALSE
  }
  test_result <- compare_means(formula, data, method = method, paired = paired, p.adjust.method = "bonferroni", var.equal = var.equal)
}


# Train a random forest or XGBoost classifier using the supplied tuning parameters. ----
train_model <- function(data, features, target, type = model_type, tune_params = NULL, folds = folds, seed = 123, return_shap=FALSE, return_cv = FALSE) {
  set.seed(seed)
  if (is.null(tune_params)) {
    tune_params <- list(
      objective = "multi:softprob",
      num_class = nlevels(data[[target]]),
      eval_metric = "mlogloss",
      max_depth = 6,
      eta = 0.1,
      gamma = 0.1
    )
  }
  # Compute class weights to account for class imbalance.
  class_weights <- data %>%
    group_by(!!sym(target)) %>%
    summarise(n = n()) %>%
    mutate(weight = max(n)/n)
  
  if (!is.null(tune_params$weight_power)) {
    class_weights <- class_weights %>%mutate(weight = weight^tune_params$weight_power)
  } else {
    class_weights$weight <- 1
  }
  if (type == "randomforest") {
    # Configure a random forest classifier.
    sampsize <- min(class_weights$n)
    randomForest(
      as.formula(paste(target, "~", paste(features, collapse = "+"))),
      data = data,
      ntree = 200,
      mtry = floor(sqrt(length(features))),
      sampsize = rep(sampsize, nlevels(data[[target]])),
      strata = data[[target]],
      classwt = setNames(class_weights$weight, levels(data[[target]])))
  } else {
    # Configure an XGBoost classifier.
    dtrain <- xgb.DMatrix(
      data = data.matrix(data[, features]),
      label = as.numeric(data[[target]]) - 1,
      weight = class_weights$weight[match(data[[target]], class_weights[[target]])]
    )
    # print(tune_params)
    # print(class_weights$weight)
    # Perform cross-validation for hyperparameter tuning.
    cv <- xgb.cv(
      params = tune_params[!names(tune_params) %in% c("nrounds", "weight_power")],
      data = dtrain,
      nrounds = tune_params$nrounds,
      nfold = 5,
      folds = folds,
      early_stopping_rounds = 10,
      verbose = 0,
      prediction = TRUE
    )
    if (return_cv) {
      return(cv)
    }
    else {
      # print(paste("best iteration:", cv$best_iteration))
      model <- xgb.train(
        params = tune_params[!names(tune_params) %in% c("nrounds", "weight_power")],
        data = dtrain,
        nrounds = cv$best_iteration,
        verbose = 0
      )
      if (return_shap) {
        shape_values <- shapviz(model, X_pred = data.matrix(data[, features]))
        return(shape_values)
      }
      return(model)
    }
  }
}


# Return predicted class probabilities and class labels using a unified interface. ----
predict_unified <- function(model, newdata, features, target_levels, type = model_type) {
  if (type == "randomforest") {
    prob <- predict(model, newdata = newdata[, features], type = "prob")
    pred <- predict(model, newdata = newdata[, features])
  } else {
    prob_raw <- predict(model, xgb.DMatrix(data.matrix(newdata[, features])))
    prob <- matrix(prob_raw, ncol = length(target_levels), byrow = TRUE)
    pred <- target_levels[max.col(prob)]
  }
  
  list(
    prob = as.data.frame(prob) %>% setNames(paste0("Prob_", target_levels)),
    pred = tibble(Pred = factor(pred, levels = target_levels))
  )
}



# Patient-level stratified data partition functions ----

# Generate diagnosis-stratified cross-validation indices without splitting an examination across folds.
split_patients_into_k_folds <- function(df, patient_col = "PatientID",
                                        diagnosis_col = "Clinical.diagnosis",
                                        k = 5, seed = 123, returnTrain = FALSE) {
  set.seed(seed)
  unique_patients <- df %>%
    distinct(across(all_of(c(patient_col, diagnosis_col))))
  if (anyDuplicated(unique_patients[[patient_col]])) {
    stop("Each examination must have a single diagnosis label.")
  }
  patient_folds <- caret::createFolds(unique_patients[[diagnosis_col]],
                                     k = k, list = TRUE, returnTrain = FALSE)
  fold_indices <- lapply(patient_folds, function(fold_patients) {
    patients_in_fold <- unique_patients[[patient_col]][fold_patients]
    which(df[[patient_col]] %in% patients_in_fold)
  })
  if (returnTrain) {
    fold_indices <- lapply(fold_indices, function(x) setdiff(seq_len(nrow(df)), x))
  }
  fold_indices
}

# Split examinations into diagnosis-stratified training and held-out testing sets.
split_preserving_diagnosis_caret <- function(df, patient_col = "PatientID",
                                             diagnosis_col = "Clinical.diagnosis",
                                             ratio = 0.7, seed = 123) {
  set.seed(seed)
  unique_patients <- df %>%
    distinct(across(all_of(c(patient_col, diagnosis_col))))
  if (anyDuplicated(unique_patients[[patient_col]])) {
    stop("Each examination must have a single diagnosis label.")
  }
  train_indices <- caret::createDataPartition(
    unique_patients[[diagnosis_col]], p = ratio, list = FALSE
  )
  train_patients <- unique_patients[[patient_col]][as.vector(train_indices)]
  df_train <- df[df[[patient_col]] %in% train_patients, ]
  df_test <- df[!df[[patient_col]] %in% train_patients, ]
  list(training = df_train, validation = df_test)
}
