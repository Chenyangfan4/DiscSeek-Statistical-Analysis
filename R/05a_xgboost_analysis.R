# DiscSeek | Exploratory XGBoost disease-subtype classification
# Outputs: Table 3, Figure 3a–c, and Supplementary Table S10 hyperparameter settings.
# The final optional block requires sensitivity-analysis outputs from notebook 05b.

source("R/00_config.R")
library(tidyverse)
library(xgboost)
library(caret)
library(multiROC)
library(shapviz)
library(patchwork)
library(boot)
library(RColorBrewer)
library(pROC)

# 1. Load shared patient-level train/test partitions ----
modeling_context <- readr::read_rds(file.path(INTERMEDIATE_DIR, "modeling_context.rds"))
list2env(modeling_context, envir = environment())

# 2. Create five patient-level cross-validation folds ----
features <- setdiff(names(patient_data), c("ScanGroup", diag_var, patient_var))  # Select model predictors from patient-level anatomical features.
folds_for_groups <- list()
for (scan_group in patient_scan_groups) {
  df <- train_val_list$balanced_training %>% filter(str_detect(ScanGroup, scan_group))
  folds_for_groups[[scan_group]] <- split_patients_into_k_folds(df, "ImageID", diag_var)  # 
  print(sapply(folds_for_groups[[scan_group]], length))
}

# 3. Table S10 | Staged XGBoost hyperparameter optimization ----
# Repeat the staged search, carrying forward the selected parameters between stages.
other_param_setting <- NULL
# Choose the hyperparameters searched at the current stage; retain earlier choices. ----
fine_tune_params <- expand.grid(
  # Stage 1: candidate learning rates.
  # eta = seq(0.01, 0.3, by=0.01),  # c(0.01, 0.05, 0.1, 0.3),
  # Stage 2: candidate tree depth and minimum child weight.
  # max_depth = c(3, 6, 9, 12),  #c(6, 9, 12),
  # min_child_weight = seq(1, 10, by=2),  
  # Stage 3: candidate splitting penalties.
  # gamma = seq(0, 5, by=0.1),  # c(0, 0.1, 0.2),
  # Stage 4: candidate column and row sampling fractions.
  # colsample_bytree = seq(0.5, 1.0, by=0.1),
  # subsample = seq(0.5, 1.0, by=0.1),
  # Stage 5: candidate regularization strengths.
  lambda = c(1e-5, 1e-2, 0.1, 1, 2, 4),
  alpha = c(1e-5, 1e-2, 0, 0.1, 1, 10),
  # Parameters held fixed during the current stage.
  nrounds = 200,
  objective = "multi:softprob",
  num_class = nlevels(patient_data[[diag_var]]),
  eval_metric = c("mlogloss")
)
cha_param <- c("objective", "eval_metric")
# cha_param <- c(cha_param, c("device", "tree_method", "sampling_method"))
# cha_param <- c("sampling_method")
fine_tune_params <- fine_tune_params %>% 
  mutate(across(cha_param, as.character))
total_params <- nrow(fine_tune_params)

# Search the selected hyperparameter grid and evaluate classification metrics. ----
progress_width <- 50  # Console progress bar width.
performance_report_patient <- map_df(patient_scan_groups, ~{
  current_group <- .x
  df <- train_val_list$balanced_training %>% filter(str_detect(ScanGroup, current_group))
  test_df <- train_val_list$balanced_validation %>% filter(str_detect(ScanGroup, current_group))
  
  cat(strrep("=", progress_width * 2), '\n')
  cat("Start process:", current_group, "\n")
  
#     y_true <- model.matrix(as.formula(paste("~", diag_var, "- 1")), data = df) %>% as.data.frame()
#     names(y_true) <- str_replace(names(y_true), diag_var, "") %>% paste0("_true")
  y_true_test <- model.matrix(as.formula(paste("~", diag_var, "- 1")), data = test_df) %>% as.data.frame()
  names(y_true_test) <- str_replace(names(y_true_test), diag_var, "") %>% paste0("_true")
  
  target_levels <- levels(df[[diag_var]])
  best_metrics <- list(
    params = NULL,
    accuracy = -1,
    balance_accuracy = -1,
    kappa = -1,
    micro_f1 = -1,
    macro_f1 = -1,
    micro_auc = -1,
    macro_auc = -1
  )
  folds <- folds_for_groups[[current_group]]
  start_time <- Sys.time()
  
  for (i in seq_len(nrow(fine_tune_params))) {
    tune_params <- fine_tune_params[i,]
    
    if (!(is.null(other_param_setting))){
      # Include hyperparameters selected in the preceding search stage.
      tune_params <- other_param_setting %>%
        filter(ScanGroup==current_group) %>%
        dplyr::select(-c(ScanGroup, intersect(names(tune_params), names(other_param_setting)))) %>%
        bind_cols(tune_params)
    }
    
    # Convert the tuning record to a named list.
    tune_params <- lapply(names(tune_params), function(x) {
      tune_params[[x]]
    }) %>% setNames(names(tune_params) %>% str_remove("param_"))
    
    # Train using five-fold cross-validation.
    model <- train_model(df, features, diag_var, model_type, 
                         tune_params = tune_params, folds = folds)
    
    # Predict on the held-out testing partition.
    preds <- predict_unified(model, test_df, features, levels(patient_data[[diag_var]]), model_type)
    names(preds$prob) <- names(preds$prob) %>% str_remove("Prob_") %>% paste0("_pred_xgb")
    
    # Evaluate classification performance.
    conf_res <- confusionMatrix(preds$pred$Pred, test_df[[diag_var]], mode = "everything")    
    TP <- conf_res$table %>% diag() %>% sum()
    FP <- sum(colSums(conf_res$table) - diag(conf_res$table))
    FN <- sum(rowSums(conf_res$table) - diag(conf_res$table)) 
    metrics <- list(
      accuracy = conf_res$overall["Accuracy"],
      kappa = conf_res$overall["Kappa"],
      balance_accuracy = mean(conf_res$byClass[,"Balanced Accuracy"]),
      micro_f1 = 2 * TP / (2 * TP + FN + FP),
      macro_f1 = mean(conf_res$byClass[, "F1"])
    )
    
    roc_res <- multi_roc(bind_cols(y_true_test, preds$prob))
    metrics$micro_auc <- roc_res$AUC$xgb$micro
    metrics$macro_auc <- roc_res$AUC$xgb$macro
    
    if (best_metrics$macro_auc < metrics$macro_auc) {
      # Select candidate parameters based on macro-averaged AUC.
      best_metrics <- modifyList(best_metrics, list(
        params = tune_params,
        accuracy = metrics$accuracy,
        kappa = metrics$kappa,
        balance_accuracy = metrics$balance_accuracy,
        micro_f1 = metrics$micro_f1,
        macro_f1 = metrics$macro_f1,
        micro_auc = metrics$micro_auc,
        macro_auc = metrics$macro_auc
      ))  
    }
    
    cat("\r[")
    progress <- round(i/total_params*progress_width)
    cat(strrep(">", progress))
    if(progress < progress_width) cat(strrep(" ", progress_width - progress))
    cat(sprintf("] %d/%d (%.1f%%) | Best Macro AUC: %.4f", 
                i, total_params, i/total_params*100, best_metrics$macro_auc))
    flush.console()
  }

  elapsed_time <- difftime(Sys.time(), start_time, units = "secs")

  cat(sprintf("\nFinish process Group %s, best Macro AUC=%.4f (cost %.1fs) \n",
              current_group, best_metrics$macro_auc, elapsed_time))
  cat(strrep("=", progress_width * 2), "\n\n")
  # print(paste("training and testing for",current_group, ":", best_metrics$macro_auc))
  #   print(names(best_metrics$params))
  result <- bind_cols(
    # Collect the selected parameter and performance summary.
    tibble(
      ScanGroup = current_group,
      Accuracy = best_metrics$accuracy,
      Kappa = best_metrics$kappa,
      Balance_Accuracy = best_metrics$balance_accuracy,
      Micro_F1 = best_metrics$micro_f1,
      Macro_F1 = best_metrics$macro_f1,
      Micro_AUC = best_metrics$micro_auc,
      Macro_AUC = best_metrics$macro_auc
    ),
    # Prefix names for selected hyperparameter columns.
    as_tibble(best_metrics$params) %>% 
      rename_with(~ paste0("param_", .))
  )
  
  return(result)
})

# Use the selected parameters as the starting point for the next stage. ----
other_param_setting <- performance_report_patient %>% 
  dplyr::select(c("ScanGroup", names(performance_report_patient)[str_detect(names(performance_report_patient), "param_")]))
names(other_param_setting) <- gsub("^param_", "", names(other_param_setting))

# 4. Table 3 | Final XGBoost models and held-out metrics ----
run_group_models <- function(return_type = c("model", "metrics", "predictions")) {
  return_type <- match.arg(return_type)
  full_result <- map(patient_scan_groups, ~{
    current_group <- .x
    # Retrieve the group-specific training and held-out testing partitions.
    df <- train_val_list$balanced_training %>% filter(str_detect(ScanGroup, current_group))
    test_df <- train_val_list$balanced_validation %>% filter(str_detect(ScanGroup, current_group))
    #   print(table(df[, diag_var]))
    # Build one-hot clinical diagnosis targets for ROC analysis.
    y_true <- model.matrix(as.formula(paste("~", diag_var, "- 1")), data = df) %>% as.data.frame()
    names(y_true) <- str_replace(names(y_true), diag_var, "") %>% paste0("_true")
    y_true_test <- model.matrix(as.formula(paste("~", diag_var, "- 1")), data = test_df) %>% as.data.frame()
    names(y_true_test) <- str_replace(names(y_true_test), diag_var, "") %>% paste0("_true")
    
    target_levels <- levels(df[[diag_var]])  # Use the observed diagnosis levels as the prediction classes.
    folds <- folds_for_groups[[current_group]]  # Load the group-specific patient-level cross-validation folds.
    
    # Retrieve hyperparameters selected in the preceding tuning stage.
    tune_params <- performance_report_patient %>% 
      filter(ScanGroup == current_group) %>% 
      dplyr::select(names(performance_report_patient)[str_detect(names(performance_report_patient), "param_")])
    
    # Convert selected parameter fields into a named list.
    tune_params <- lapply(names(tune_params), function(x) {
      tune_params[[x]]
    }) %>% setNames(names(tune_params) %>% str_remove("param_"))
    #   print(as.data.frame(tune_params)) 
    
    # Train the classifier on the resampled training partition.
    model <- train_model(df, features, diag_var, model_type, 
                         tune_params = tune_params, folds = folds)
    
    switch(return_type,
           "model" = {list(
             model = model, 
             scan_group = current_group
           )
           },
           "metrics" = {
             calculate_metrics <- function(model, data, true_values, require_ci = FALSE) {
               preds <- predict_unified(model, data, features, target_levels, model_type)
               conf_res <- confusionMatrix(preds$pred$Pred, data[[diag_var]], mode = "everything")
               
               TP <- diag(conf_res$table) %>% sum()
               FP <- sum(colSums(conf_res$table) - diag(conf_res$table))
               FN <- sum(rowSums(conf_res$table) - diag(conf_res$table))
               
               names(preds$prob) <- names(preds$prob) %>% str_remove("Prob_") %>% paste0("_pred_xgb")
               roc_res <- multi_roc(bind_cols(true_values, preds$prob))
               
               results <- list()
               results$metrics <- list(
                 accuracy = unname(conf_res$overall["Accuracy"]),
                 kappa = unname(conf_res$overall["Kappa"]),
                 balance_accuracy = mean(conf_res$byClass[,"Balanced Accuracy"] %>% replace(., is.na(.), 0), na.rm = TRUE),
                 macro_f1 = mean(conf_res$byClass[,"F1"] %>% replace(., is.na(.), 0), na.rm = TRUE),
                 micro_f1 = 2 * TP / (2 * TP + FN + FP),
                 micro_auc = roc_res$AUC$xgb$micro,
                 macro_auc = roc_res$AUC$xgb$macro
               )
               if (!require_ci) {
                 return(results$metrics)
               } else {
                 # Bootstrap accuracy, balanced accuracy, and F1 confidence intervals.
                 calculate_acc_res <- function(input, indice) {
                   input <- input[indice,]
                   conf_res <- confusionMatrix(input$Pred, input[[diag_var]], mode = "everything")
                   TP <- diag(conf_res$table) %>% sum()
                   FP <- sum(colSums(conf_res$table) - diag(conf_res$table))
                   FN <- sum(rowSums(conf_res$table) - diag(conf_res$table))
                   return(c(
                     accuracy = unname(conf_res$overall["Accuracy"]),
                     balance_accuracy = mean(conf_res$byClass[,"Balanced Accuracy"] %>% replace(., is.na(.), 0), na.rm = TRUE),
                     macro_f1 = mean(conf_res$byClass[,"F1"] %>% replace(., is.na(.), 0), na.rm = TRUE),
                     micro_f1 = 2 * TP / (2 * TP + FN + FP)
                   ))
                 }
                 ci_res <- boot(data = bind_cols(preds$pred, data[,diag_var]), statistic = calculate_acc_res, R = 1000, strata = data[[diag_var]])
                 results$ci_lower <- list(
                   accuracy = boot.ci(ci_res, type = "perc", index = 1)$percent[1,4] ,
                   kappa = NA,
                   balance_accuracy = boot.ci(ci_res, type = "perc", index = 2)$percent[1,4],
                   macro_f1 = boot.ci(ci_res, type = "perc", index = 3)$percent[1,4],
                   micro_f1 = boot.ci(ci_res, type = "perc", index = 4)$percent[1,4])
                 results$ci_upper <- list(
                   accuracy = boot.ci(ci_res, type = "perc", index = 1)$percent[1,5],
                   kappa = NA,
                   balance_accuracy = boot.ci(ci_res, type = "perc", index = 2)$percent[1,5],
                   macro_f1 = boot.ci(ci_res, type = "perc", index = 3)$percent[1,5],
                   micro_f1 = boot.ci(ci_res, type = "perc", index = 4)$percent[1,5]
                 )
                 print(names(results$ci_lower))
                 # Bootstrap micro- and macro-averaged AUC confidence intervals.
                 calculate_auc_res <- function(input, indice) {
                   roc_res <- multi_roc(input[indice,])
                   return(c(
                     micro_auc = roc_res$AUC$xgb$micro,
                     macro_auc = roc_res$AUC$xgb$macro
                   ))
                 }
                 ci_res <- boot(data = bind_cols(true_values, preds$prob), statistic = calculate_auc_res, R = 1000, strata = data[[diag_var]])
                 results$ci_lower <- c(results$ci_lower, list(
                   micro_auc = boot.ci(ci_res, type = "perc", index = 1)$percent[1,4],
                   macro_auc = boot.ci(ci_res, type = "perc", index = 2)$percent[1,4]
                 ))
                 results$ci_upper <- c(results$ci_upper, list(
                   micro_auc = boot.ci(ci_res, type = "perc", index = 1)$percent[1,5],
                   macro_auc = boot.ci(ci_res, type = "perc", index = 2)$percent[1,5]
                 ))
                 
                 return(results)
               }
             }
             # train_metrics <- calculate_metrics(model, df, y_true)
             test_metrics <- calculate_metrics(model, test_df, y_true_test, require_ci=TRUE)
             
             # Optional: format test metrics as display strings with confidence intervals.
             # test_metrics <- map_dfr(test_metrics, as.data.frame) %>% 
             #   mutate(ID = names(test_metrics)) %>% 
             #   column_to_rownames(var = "ID") %>% 
             #   t() %>% as.data.frame() %>% 
             #   mutate(formatted = sprintf("%.3f (%.2f, %.2f)", metrics, ci_lower, ci_upper)) %>%
             #   dplyr::select(formatted) %>%
             #   t() %>% as.data.frame() %>% remove_rownames()
             
             # Keep numerical confidence interval bounds for downstream analysis.
             test_metrics <- map_dfr(test_metrics, as.data.frame) %>% 
               mutate(ID = names(test_metrics)) %>% 
               column_to_rownames(var = "ID") %>% 
               t() %>% as.data.frame() %>% 
               mutate(ScanGroup = current_group, Dataset = "test") %>%
               rownames_to_column() %>%
               rename(Metrics = rowname, Value = metrics)
             
             # bind_cols(tibble(ScanGroup = current_group),
             #           as_tibble(test_metrics) %>%
             #             rename_with(~ paste0("test_", .)),
             #           as_tibble(train_metrics) %>%
             #             rename_with(~ paste0("train_", .)))
           },
           "predictions" = {
             preds <- predict_unified(model, df, features, levels(patient_data[[diag_var]]), model_type)
             train_preds <- preds$prob %>% 
               mutate(Dataset = "train") %>% 
               bind_cols(df[,diag_var])
             
             preds <- predict_unified(model, test_df, features, levels(patient_data[[diag_var]]), model_type)
             test_preds <- preds$prob %>% 
               mutate(Dataset = "test") %>% 
               bind_cols(test_df[,diag_var])
             # Return class-probability predictions.
             bind_rows(train_preds, test_preds) %>% mutate(ScanGroup = current_group)
           }
    )
  })

return(full_result)
}

target_levels <- levels(patient_data[[diag_var]])
# Return final trained models for each planning group.
full_result <- run_group_models(return_type = "model")
full_result <- setNames(full_result, patient_scan_groups)
# Calculate held-out classification metrics.
formated_performance_report <- run_group_models(return_type = "metrics")  
formated_performance_report <- bind_rows(formated_performance_report)
# Generate training and held-out prediction probabilities.
total_preds <- run_group_models(return_type = "predictions")  
total_preds <- bind_rows(total_preds)
total_preds$Pred <- target_levels[max.col(total_preds[,1:3])]
total_preds$Pred <- factor(total_preds$Pred, levels = target_levels)

# Export selected hyperparameters and test-set metrics.
write_csv(performance_report_patient, file.path(MODELING_TABLE_DIR, paste0(model_type, "_patient_metrics_split(smote).csv")))
# Export held-out classification metrics and confidence intervals.
write_csv(formated_performance_report, file.path(MODELING_TABLE_DIR, paste0(model_type, "_patient_metrics(formate)_split(smote).csv")))

# 5. Figure 3a | Held-out test ROC curves ----
# Compute class-wise and averaged ROC curve coordinates. ----
total_roc_df <- total_preds %>% distinct(ScanGroup) %>% pmap_dfr(~ {
  current_group <- ..1
  df <- total_preds %>% filter(ScanGroup == current_group, Dataset == "test")
  conf_res <- confusionMatrix(df$Pred, df[[diag_var]], mode = "everything")
  
  prob <- df[,1:3]
  names(prob) <- names(prob) %>% str_remove("Prob_") %>% paste0("_pred_xgb")
  y_true <- model.matrix(as.formula(paste("~", diag_var, "- 1")), data = df) %>% as.data.frame()
  names(y_true) <- str_replace(names(y_true), diag_var, "") %>% paste0("_true")
  roc_res <- multi_roc(bind_cols(y_true, prob))
  roc_res_df <- NULL
  for (cls in names(roc_res$Specificity$xgb)){
    roc_res_df <- roc_res_df %>% bind_rows(
      as.data.frame(
        list(
          Sensitivity = roc_res$Sensitivity$xgb[[cls]], 
          Specificity = roc_res$Specificity$xgb[[cls]],
          Type = cls
        )
      )
    )
  }
  roc_res_df <- roc_res_df %>% mutate(Type = factor(Type), ScanGroup = current_group)
})

# Plot test-set ROC curves. ----
thresh_df <- as.data.frame(list(
  x = 0,
  y = 0,
  x_end = 1,
  y_end = 1,
  Type = unique(total_roc_df$Type) %>% 
    str_replace_all(., c("^micro$" = "Micro",
                         "^macro$" = "Macro"))
))
label_df <- as.data.frame(list(
  x = 1.0,
  y = 0.025,
  Type = unique(total_roc_df$Type) %>% 
    str_replace_all(., c("^micro$" = "Micro",
                         "^macro$" = "Macro"))
)) %>% mutate(label = paste(Type, "ROC"))
roc_plot <- total_roc_df %>% 
  mutate(Type = str_replace_all(Type, c("^micro$" = "Micro", "^macro$" = "Macro")),
         ScanGroup = factor(ScanGroup, levels = patient_scan_groups[c(4,3,1,2)], ordered = TRUE)) %>%
  mutate(ScanGroup = fct_relabel(ScanGroup, ~{
    str_replace_all(., c(
      "^Manual$" = "Multi-center Manual",
      "^DT_UIH_AI$" = "H1/S1 DiscSeek",
      "^DT_UIH_Interactive$" = "H1/S1 DiscSeek-assisted",
      "^DT_UIH_Manual$" = "H1/S1 Manual"
    ))})) %>%
  ggplot(aes(y = Sensitivity, x = 1-Specificity, group = ScanGroup, color = ScanGroup)) + 
  geom_line() + 
  # scale_color_brewer(palette = "Set2", type = "qual") + 
  scale_color_manual(values = color_palettes) + 
  # scale_linetype_manual(values = setdiff(1:6,3)) + 
  geom_segment(aes(x = x, xend = x_end, y = y, yend = y_end), 
               data = thresh_df , inherit.aes = FALSE, colour = "grey50", linetype = 2) + 
  geom_label(aes(x = x, y = y, label = label), data = label_df, 
             size = 8, size.unit = "pt",
             inherit.aes = FALSE, hjust = 1, vjust = 0) + 
  facet_wrap(~ Type, nrow = 1, dir = "v", 
             # labeller = labeller(Type = function(value) paste0("ROC(", value, ")"))
  ) + 
  labs(
    # title = "Receiver Operating Characteristic", 
    x = "False Positive Rate", y = "True Positive Rate",
    color = "Planning Group", linetype = NULL) + 
  theme(plot.title = element_blank(),
        # legend.position = c(0.5, 1.1),
        aspect.ratio = 1.0,
        strip.text = element_blank(),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.25),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.25),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
        axis.title = element_text(size = 10),
        axis.text = element_text(color = "black", size = 6),
        legend.box.spacing = unit(0,'pt'),
        legend.text = element_text(size = 8),
        legend.position = c(0.5, 1.1),  # "top",
        legend.direction = "horizontal",
        legend.key = element_rect(fill = "white", colour = "grey", linewidth = 0),
        # legend.title.position = "left",
        # legend.justification = "center",
        # legend.margin = ggplot2::margin(0,0,0,0),
        legend.box.margin = ggplot2::margin(0,0,0,0),
        legend.key.height = unit(0.1, "mm"),
        # legend.title = element_text(margin = ggplot2::margin(r=10)),
        plot.margin = ggplot2::margin(0,0,0,0, unit = "pt")
  )
ggsave(file.path(MODELING_FIGURE_DIR, "xgboost_roc.png"), roc_plot, width = 7, height = 2.2)

# 6. Figure 3b | Classification performance and AUC comparisons ----
# Compare micro-averaged AUC between planning groups.
calc_micro_auc_roc <- function(df, target_levels) {
  # Construct one-versus-rest multiclass data for ROC testing.
  y_true <- unlist(lapply(target_levels, function(cls) {
    ifelse(df$Clinical.diagnosis == cls, 1, 0)
  }))
  
  y_pred <- unlist(lapply(target_levels, function(cls) {
    df[[paste0("Prob_", cls)]]
  }))
  
  # Calculate ROC curves for paired comparisons where applicable.
  pROC::roc(response = y_true, predictor = y_pred)
}

total_rocs <- total_preds %>% 
  distinct(ScanGroup, Clinical.diagnosis) %>% 
  pmap(~ {
    current_group <- ..1
    diag_type <- ..2
  y_true <- total_preds %>% 
    filter(Dataset == "test", ScanGroup == current_group) %>% 
    mutate(label = ifelse(!!sym(diag_var) == diag_type, 1, 0),
           pred = !!sym(paste0("Prob_", diag_type))) %>%
    dplyr::select(c(label, pred))
  list(
    roc = pROC::roc(response = y_true[, "label"], predictor = y_true[,"pred"]),
    Clinical.diagnosis = as.character(diag_type),
    ScanGroup = as.character(current_group)
  )
}) %>%
  set_names(map_chr(., ~ paste(.x$ScanGroup, .x$Clinical.diagnosis, sep = "&")))

for (current_group in unique(total_preds$ScanGroup)) {
  total_rocs[[paste(current_group, "Micro", sep = "&")]] <-   list(
    roc = calc_micro_auc_roc(total_preds %>% 
                               filter(Dataset == "test", ScanGroup == current_group),
                             target_levels),
    Clinical.diagnosis = "Micro",
    ScanGroup = as.character(current_group)
  )
}

group_paires <- combn(patient_scan_groups, 2, simplify = FALSE)

pairwise_aus_comp <- pmap_dfr(list(seq(1, length(group_paires))), function(i) {
  group1 <- group_paires[[i]][1]
  group2 <- group_paires[[i]][2]
  print(paste(group1, group2, sep = ", "))
  res_results <- map_dfr(c(target_levels, "Micro"), function(diag_type) {
    res <- pROC::roc.test(
      roc1 = total_rocs[[paste0(group1, "&", diag_type)]]$roc,
      roc2 = total_rocs[[paste0(group2, "&", diag_type)]]$roc,
      method = "delong",
      paired = FALSE
    )
    tibble(
      Group1 = group1,
      Group2 = group2,
      Diagnosis = diag_type,
      Method = "delong",
      statistic = res$statistic,
      p_value = res$p.value,
      AUC1 = res$estimate["AUC of roc1"],
      AUC2 = res$estimate["AUC of roc2"]
    )
  })
})

# Estimate macro-AUC differences using bootstrap resampling.
calc_macro_auc <- function(input, indices) {
  if (!missing(indices)) {
    input <- input[indices, ]
  }
  roc_res <- multi_roc(input)
  return(roc_res$AUC$xgb$macro)
}

n_bootstraps <- 1000

pairwise_aus_comp <- pmap_dfr(list(seq(1, length(group_paires))), function(i) {
  group1 <- group_paires[[i]][1]
  group2 <- group_paires[[i]][2]
  print(paste(group1, group2, sep = ", "))
  prob_1 <- total_preds %>% 
    filter(ScanGroup == group1, Dataset == "test") %>% 
    dplyr::select(paste0("Prob_", target_levels))
  names(prob_1) <- names(prob_1) %>% str_remove("Prob_") %>% paste0("_pred_xgb")
  prob_2 <- total_preds %>% 
    filter(ScanGroup == group2, Dataset == "test") %>% 
    dplyr::select(paste0("Prob_", target_levels))
  names(prob_2) <- names(prob_2) %>% str_remove("Prob_") %>% paste0("_pred_xgb")
  
  y_true_1 <- model.matrix(as.formula(paste("~", diag_var, "- 1")), data =  total_preds %>% 
                                filter(ScanGroup == group1, Dataset == "test")) %>% as.data.frame()
  names(y_true_1) <- str_replace(names(y_true_1), diag_var, "") %>% paste0("_true")
  
  y_true_2 <- model.matrix(as.formula(paste("~", diag_var, "- 1")), data =  total_preds %>% 
                             filter(ScanGroup == group2, Dataset == "test")) %>% as.data.frame()
  names(y_true_2) <- str_replace(names(y_true_2), diag_var, "") %>% paste0("_true")
  
  df_filter_1 <- bind_cols(prob_1, y_true_1)
  df_filter_2 <- bind_cols(prob_2, y_true_2)
  
  
  auc1_obs <- calc_macro_auc(df_filter_1)
  auc2_obs <- calc_macro_auc(df_filter_2)
  delta_auc_obs <- auc1_obs - auc2_obs
  
  boot_result1 <- boot(data = df_filter_1, statistic = calc_macro_auc, R = n_bootstraps, 
                       strata = total_preds %>% 
                         filter(ScanGroup == group1, Dataset == "test") %>% 
                         pull(!!sym(diag_var)))
  boot_result2 <- boot(data = df_filter_2, statistic = calc_macro_auc, R = n_bootstraps, 
                       strata = total_preds %>% 
                         filter(ScanGroup == group2, Dataset == "test") %>% 
                         pull(!!sym(diag_var)))
  delta_aucs <- boot_result1$t - boot_result2$t
  
  ci_result <- boot.ci(boot.out = list(t = delta_aucs, R = n_bootstraps, t0=delta_auc_obs), 
                       type = "perc")  
  
  prop_greater <- mean(delta_aucs > 0, na.rm = TRUE)
  prop_less <- mean(delta_aucs < 0, na.rm = TRUE)
  p_value <- 2 * min(prop_greater, prop_less)
  tibble(
    Group1 = group1,
    Group2 = group2,
    Diagnosis = "Macro",
    Method = "Bootstrap",
    AUC1 = auc1_obs,
    AUC2 = auc2_obs,
    p_value = p_value,
    statistic = NA
  )
}) %>% bind_rows(pairwise_aus_comp)

write_csv(pairwise_aus_comp, file.path(MODELING_TABLE_DIR, "pairwise_comparision_auc.csv"))

# Generate the classification metric plots.
colors <- c("#3951A2", "#72AACF", "#CAE8F2", "#FEFBBA", "#FDB96B", "#EC5D3B", "#A80326")
color_palettes <- c("#FFB117", "#8979D0","#6DB50A", "#5F9BED")
color_palettes <- c("#FFA725","#FF6B62","#BE588D","#444577")
color_palettes <- c("#40004B", "#8B559B", "#CAAFD4", "#C2E6BC", "#56AA5E", "#095624")[c(2,3,4,5)]
color_palettes <- c(colors)[c(1,2,5,6)]
visible_metrics <- c("balance_accuracy", "macro_f1", "micro_f1", "macro_auc", "micro_auc")
# Plot test-set classification metrics.
stat_df <- formated_performance_report %>% 
  # pivot_longer(
  #   cols = -ScanGroup,
  #   names_to = c("Dataset", "Metrics"),
  #   values_to = "Value"
  # ) %>% 
  filter(Metrics %in% visible_metrics, Dataset=="test") %>%
  mutate(Metrics = factor(Metrics, levels = visible_metrics, ordered = TRUE), 
         ScanGroup = factor(ScanGroup, levels = patient_scan_groups[c(4,3,1,2)], ordered = TRUE),
         # Type = "Quantativce Classification Performance"
  ) %>%
  mutate(
    Metrics = fct_relabel(Metrics, ~ {
      str_replace_all(., "_", " ") %>% 
        str_to_title() %>% 
        str_replace_all(., "Balance", "Balanced") %>%
        str_replace_all(., "Auc", "AUC")
    }),
    ScanGroup = fct_relabel(
      ScanGroup, ~ {
        str_replace_all(
          ., c(
            "^Manual$" = "Multi-center Manual",
            "^DT_UIH_AI$" = "H1/S1 DiscSeek",
            "^DT_UIH_Interactive$" = "H1/S1 DiscSeek-assisted",
            "^DT_UIH_Manual$" = "H1/S1 Manual"
          ))
      })
  )

metric_level_names <- levels(stat_df$Metrics)
bar_plot <- stat_df %>% 
  ggplot(aes(y = Value, x = Metrics, fill = ScanGroup)) + 
  geom_bar(width = 0.7, position = position_dodge(0.8), stat = "identity", show.legend = TRUE, colour = "grey40", linewidth = 0.25) + 
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                width = 0.2,
                linewidth = 0.25,
                position = position_dodge(0.8)) + 
  # facet_wrap(~ Dataset, nrow = 1, scales = "free_x", labeller = "label_both") + 
  # facet_wrap(~ Type) + 
  labs(x = "Diagnostic Metric", y = "Metric Value", fill = "Planning Group", 
       # title = "Classification Performance Comparison among Scanning Methods"
       ) + 
  # scale_fill_brewer(palette = "Set2", type = "qual") + 
  scale_fill_manual(values = color_palettes) + 
  # scale_y_continuous(breaks = seq(0,1,by=0.2)) + 
  # coord_cartesian(ylim = c(0, 1.0), xlim = c(0.59,5.41), expand = FALSE) + 
  theme_linedraw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        # aspect.ratio = 2.0,
        plot.margin = ggplot2::margin(0,0,0,0, unit = "mm"),
        axis.title = element_text(size = 10),
        axis.text = element_text(color = "black", size = 8),
        # axis.text.x = element_text(color = "black", size = 10),
        panel.grid.major.x = element_blank(),
        legend.box.spacing = unit(0,'pt'),
        legend.text = element_text(size = 8),
        legend.position = "top",
        legend.direction = "horizontal",
        # legend.title.position = "left",
        # legend.justification = "center",
        legend.margin = ggplot2::margin(0,0,0,0),
        legend.key.height = unit(8,"pt"),
        legend.box.margin = ggplot2::margin(0,0,0,0),
        # legend.title = element_text(margin = ggplot2::margin(r=10))
        )

# Determine bar coordinates for macro- and micro-AUC comparisons.
gg_data <- ggplot_build(bar_plot)$data[[1]]
x_pos <- list()
for (metric_type in levels(stat_df$Metrics)) {
  x_pos[[metric_type]] <- gg_data %>% 
    filter((x < match(metric_type, metric_level_names) + 0.5) & (x > match(metric_type, metric_level_names) - 0.5)) %>% 
    pull(x)%>%unique() %>% sort()
}
# Position pairwise AUC significance annotations.
posthoc_df <- pairwise_aus_comp %>% 
  filter(Diagnosis %in% c("Micro", "Macro")) %>%
  mutate(Metrics = paste(Diagnosis, "AUC")) %>%
  rename(group1 = Group1, group2 = Group2) %>%
  mutate(p.signif = case_when(p_value < 0.0001 ~ "****", 
                              p_value < 0.001 ~ "***", 
                              p_value < 0.01 ~ "**", 
                              p_value < 0.05 ~ "*", 
                              TRUE ~ "ns"),
         group1 = factor(group1, levels = patient_scan_groups[c(4,3,1,2)], ordered = TRUE),
         group2 = factor(group2, levels = patient_scan_groups[c(4,3,1,2)], ordered = TRUE)
         ) %>%
  mutate(
    group1 = fct_relabel(
      group1, ~ {
        str_replace_all(
          ., c(
            "^Manual$" = "Multi-center Manual",
            "^DT_UIH_AI$" = "H1/S1 DiscSeek",
            "^DT_UIH_Interactive$" = "H1/S1 DiscSeek-assisted",
            "^DT_UIH_Manual$" = "H1/S1 Manual"
          ))
      }),
    group2 = fct_relabel(
      group2, ~ {
        str_replace_all(
          ., c(
            "^Manual$" = "Multi-center Manual",
            "^DT_UIH_AI$" = "H1/S1 DiscSeek",
            "^DT_UIH_Interactive$" = "H1/S1 DiscSeek-assisted",
            "^DT_UIH_Manual$" = "H1/S1 Manual"
          ))
      })
  )

max_y <- stat_df %>% group_by(Metrics) %>% summarise(max_y = max(ci_upper), .groups = "drop")  # Determine the maximum confidence interval bound for each metric.
y_step <- max(max_y$max_y) # Define the vertical spacing for significance annotations.

posthoc_df <- posthoc_df %>% 
  left_join(max_y, by = c("Metrics")) %>%
  mutate(Metrics = factor(Metrics, levels = metric_level_names))


posthoc_df <- posthoc_df %>%
  group_by(Metrics) %>% 
  mutate(x = map(Metrics, ~ x_pos[[.]]),
         group = list(sort(unique(c(group1, group2)), 
                           method = "radix")[order(match(sort(unique(c(group1, group2)), 
                                                              method = "radix"), 
                                                         metric_level_names))]),
         id1 = match(group1, group[[1]]),  # Map compared groups to their plotted positions.
         id2 = match(group2, group[[1]]),
         x_1 = map2_dbl(x, id1, ~ .x[.y]),
         x_2 = map2_dbl(x, id2, ~ .x[.y])) %>% ungroup()

posthoc_df <- posthoc_df %>% 
  filter(p.signif != "ns") %>%   # Suppress nonsignificant comparison annotations.
  group_by(Metrics) %>%
  mutate(
    x_min = pmin(x_1, x_2),
    x_max = pmax(x_1, x_2),
    x_text = (x_min + x_max) / 2,
    y = max_y + (row_number()) * 0.1 * y_step, 
  ) %>% ungroup() %>%
  dplyr::select(-c(x, id1, id2, x_1, x_2, max_y, group)) %>%
  arrange(Metrics, group1, group2)

# Export the version without significance annotations.
ggsave(file.path(MODELING_FIGURE_DIR, "xgboost_bar.png"), bar_plot, width = 7, height = 2.2)

# Export the version with significance annotations.
sig_bar_plot <- bar_plot +
  geom_segment(aes(x = x_min, xend = x_max, y = y, yend = y), 
               data = posthoc_df , inherit.aes = FALSE, linewidth = 0.25) + 
  geom_segment(aes(x = x_min, xend = x_min, y = 0.99 * y, yend = y), 
               data = posthoc_df , inherit.aes = FALSE, linewidth = 0.25) + 
  geom_segment(aes(x = x_max, xend = x_max, y = 0.99 * y, yend = y), 
               data = posthoc_df , inherit.aes = FALSE, linewidth = 0.25) + 
  geom_text(aes(label = p.signif, x = x_text, y = y),
            data = posthoc_df, inherit.aes = FALSE, size = 3, vjust = 0, nudge_y = -0.02) + 
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.1))
  )
ggsave(file.path(MODELING_FIGURE_DIR, "xgboost_bar_sig.eps"), sig_bar_plot,width = 174, height = 180*4.5/16, units = "mm", dpi = 300)


# Combine Figure 3a ROC and Figure 3b performance panels.
combined_plot <- wrap_plots(
  list(
    "roc" = roc_plot + 
      labs(tag = "a", x = "1-Specificity", y = "Sensitivity") +
      theme(
        legend.key.height = unit(8, "pt"),
        plot.tag = element_text(size = 10, face = "bold"),
        legend.box.spacing = unit(3,"pt"),
        legend.position = c(0.5, 1.12),
        legend.margin = ggplot2::margin(0,0,0,0),
        legend.box.margin = ggplot2::margin(0,0,0,0, unit = "pt"),
        legend.title = element_text(size = 10),
        # plot.tag.position = c(0.5, 0.07)
        ),
    "bar" = sig_bar_plot + 
      labs(tag = "b") +
      theme(
        legend.key.height = unit(8, "pt"),
        plot.tag = element_text(size = 10, face = "bold"),
        legend.box.spacing = unit(3,"pt"),
        legend.box.margin = ggplot2::margin(t=2, unit = "pt"),
        legend.position = c(0.5, 1.1),
        legend.title = element_text(size = 10)
        # plot.tag.position = c(0.5, -0.03)
        )), 
  heights = c(1.0,1.0), ncol = 1)
ggsave(file.path(MODELING_FIGURE_DIR, "xgboost_accuracy_vertical_5.png"), combined_plot, width = 180, height = 180*9/16, units = "mm", dpi = 300)

# 7. Figure 3c | XGBoost feature importance ----
# Standardize anatomical feature labels.
parameters <- c(
  "Grade.of.left.foraminal.stenosis" = "Foraminal Stenosis (L)",
  "Grade.of.right.foraminal.stenosis" = "Foraminal Stenosis (R)",
  "Height.of.left.lateral.recess..mm." = "Lateral Recess Height (L)",
  "Height.of.right.lateral.recess..mm." = "Lateral Recess Height (R)",
  "Morphology.symmetry.of.bilateral.lateral.recesses" = "Lateral Recess Symmetry",
  "Grade.of.IVD.protrusion" = "Disc Protrusion",
  "Grade.of.spinal.stenosis" = "Spinal Stenosis",
  "Sagittal.diameter.of.spinal.canal..cm." = "Sagittal Canal Diameter",
  "Transverse.diameter.of.spinal.canal..cm." = "Transverse Canal Diameter",
  "Torg.ratio" = "Torg Ratio",
  "Sagittal.diameter.of.IVD..cm." = "Sagittal Disc Diameter",
  "Transverse.diameter.of.IVD..cm." = "Transverse Disc Diameter",
  "Height.of.IVD..cm." = "Disc Height",
  "Sagittal.diameter.of.VB..cm." = "Sagittal vertebral body diameter"
)
# Compute feature importance values from the fitted classifiers.
cls_spec <- TRUE
feat_import_df <- map_dfr(patient_scan_groups, ~ {
  current_group <- .x
  model <- full_result[[current_group]]$model
  if (cls_spec) {
  import_df_1 <- xgb.importance(feature_names = features, 
                                model = model,
                                trees = seq(from=0, by=length(target_levels), length.out = model$niter)) %>%
    as.data.frame() %>% 
    mutate(ScanGroup = full_result[[current_group]]$scan_group,
           Diagnosis = target_levels[1])
  import_df_2 <- xgb.importance(feature_names = features, 
                                model = model,
                                trees = seq(from=1, by=length(target_levels), length.out = model$niter)) %>%
    as.data.frame() %>% 
    mutate(ScanGroup = full_result[[current_group]]$scan_group,
           Diagnosis = target_levels[2])
  import_df_3 <- xgb.importance(feature_names = features, 
                                model = model,
                                trees = seq(from=2, by=length(target_levels), length.out = model$niter)) %>%
    as.data.frame() %>% 
    mutate(ScanGroup = full_result[[current_group]]$scan_group,
           Diagnosis = target_levels[3])
  
  bind_rows(import_df_1, import_df_2, import_df_3)
  } else {
    xgb.importance(feature_names = features, model = model) %>% 
      as.data.frame() %>%
      mutate(ScanGroup = full_result[[current_group]]$scan_group)
  }
})
feat_import_df <- feat_import_df %>% 
  filter(!Feature %in% other_params) %>% 
  separate(Feature, c("Parameter", "IVDType"), sep = "_")

# Select the feature importance measure for visualization.
importance_col <- "Gain"  # Use Gain as the default feature importance measure; Cover and Frequency are alternatives.
feat_import_df$Value <- feat_import_df[[importance_col]]
feat_import_df$Selected_status <- 'Selected'
labels <- c("Non-selected", 'Selected')  # Label features excluded from fitted models.
scan_group_dict <- setNames(c("Manual", "H1/S1 Manual", "H1/S1 DiscSeek", "H1/S1 DiscSeek-assisted"), 
                            patient_scan_groups[c(4,3,1,2)])
plot_list <- list()
for (i in c(4,3,1,2)) {
  # Restore missing features as unselected in the complete feature grid.
  group_feat_import_df <- feat_import_df %>%
    filter(Parameter %in% names(parameters), ScanGroup == patient_scan_groups[i]) %>% 
    mutate(Selected_status = factor(Selected_status, levels = labels)) %>% 
    complete(IVDType = unique(IVDType),
             Parameter = names(parameters),
             Diagnosis = target_levels,
             fill = list(Value = NA, Selected_status = "Non-selected")) %>% # labels[1])) %>%
    mutate(Parameter = factor(Parameter, levels = names(parameters), ordered = TRUE),
           IVDType = str_replace(IVDType, "([A-Z])(\\d)([A-Z])(\\d)", "\\1\\2/\\3\\4"),
           ScanGroup = factor(ScanGroup, levels = patient_scan_groups[c(4,3,1,2)], ordered = TRUE)) %>%
    mutate(ScanGroup = fct_relabel(ScanGroup, ~{
      str_replace_all(., c(
        # "^Manual$" = "Manual(multi-center)",
        "^DT_UIH_AI$" = "H1/S1 DiscSeek",
        "^DT_UIH_Interactive$" = "H1/S1 DiscSeek-assisted",
        "^DT_UIH_Manual$" = "H1/S1 Manual"
      ))}))
  # Plot feature importance across planning groups.
  plot_list[[patient_scan_groups[i]]] <- group_feat_import_df %>%
    ggplot(aes(x = IVDType, y = Parameter)) + 
    geom_tile(aes(fill = Value), color = "black", show.legend = TRUE) + 
    geom_point(  # Mark anatomical parameters not used by the fitted classifier.
      aes(shape = Selected_status), 
      data = group_feat_import_df %>% filter(is.na(Value)), 
      size = 2, color = "grey40", stroke = 0.25
    ) + 
    facet_wrap(~ Diagnosis, nrow = 1) +
    guides(
      shape = guide_legend(
        title = NULL,
        override.aes = list(
          fill = "grey99",
          shape = 4,
          color = "black",
          size = 2),
        label = "Non-selected",
        order = 2)
    ) +
    scale_fill_gradientn(
      name = "Parameter Importance",
      colors = colors,  # brewer.pal(n=9, name = "OrRd"),  #colors[4:7],
      limits = c(0, max(feat_import_df$Value, na.rm = TRUE)),
      # values = scales::rescale(quantile(feat_import_df$Value,
      #                                   prob = seq(0, 1, length.out = length(colors)),
      #                                   na.rm = TRUE)),
      breaks = c(0, max(feat_import_df$Value, na.rm = TRUE)),
      labels = c("Low", "High"),
      na.value = "grey99"
    ) +
    scale_y_discrete(labels = parameters) + 
    scale_shape_manual(values = 4, breaks = "Non-selected", labels = "Non-selected") + 
    labs(x = "Disc Level", y = "Cervical Parameter", title = scan_group_dict[patient_scan_groups[i]]) + 
    coord_cartesian(expand = FALSE) + 
    # coord_fixed() + 
    theme(
      plot.title = element_text(size = 8, hjust = 0.5, margin = ggplot2::margin(t=5, b=3, unit = "pt")),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(size = 8, margin = ggplot2::margin(t=0,b=0)),
      legend.title = element_text(size = 10, margin = ggplot2::margin(r=10, unit = "pt")),  # , angle = 90, hjust = 0.5),
      legend.title.position = "left",
      legend.position = "bottom",
      legend.key.height = unit(5, "pt"),
      legend.box.spacing = unit(3, "pt"),
      legend.margin = ggplot2::margin(0,0,0,0),
      legend.justification = "center",
      legend.box.just = "center",
      axis.text = element_text(color = "black", size = 8),
      axis.text.x = element_text(angle = 30, hjust = 1, size = 6),
      axis.title = element_text(size = 10),
      panel.grid.minor = element_line(colour = "grey", linewidth=0.25),
      panel.grid.major = element_line(colour = "grey", linewidth = 0.25),
      panel.background = element_rect(fill = "white", colour = "grey", linewidth = 0.5),
    )
}

# Combine and export the feature importance heatmaps.
ggsave(
  file.path(
    FIGURE_DIR, 
    "modeling_prediction", 
    "xgboost_feature_importance.eps"), 
  wrap_plots(plot_list, nrow = 2, guides = "collect", 
             axes = "collect", axis_titles = "collect") &
    # plot_annotation(tag_levels = "a") & 
    theme(
      plot.tag = element_text(size = 10, margin = ggplot2::margin(0,0,0,0)),
      plot.tag.location = "plot",
      plot.margin = ggplot2::margin(l=3,r=3,unit = "pt"),
      legend.title = element_text(size = 10, margin = ggplot2::margin(r=10, unit = "pt")),  # , angle = 90, hjust = 0.5),
      legend.title.position = "left",
      legend.position = "bottom",
      legend.key.height = unit(5, "pt"),
      legend.box.spacing = unit(3, "pt"),
      legend.margin = ggplot2::margin(0,0,0,0),
      legend.justification = "center",
      legend.box.just = "center"
      ), 
  width = 180, height = 120, units = "mm", dpi = 300
  )

# Identify class-specific top features for optional manual figure annotation.
# Optional: annotate the exported feature importance panel in an external graphics editor.
top_feat_import_df <- feat_import_df %>% 
  group_by(ScanGroup, Diagnosis) %>% 
  arrange(desc(Value), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup()

# 8. Optional | Visualize sensitivity-classifier results after notebook 05b ----
# Set to TRUE after running the 05b notebook to enable this additional comparison.
RUN_OPTIONAL_CLASSIFIER_COMPARISON <- FALSE
if (RUN_OPTIONAL_CLASSIFIER_COMPARISON) {
# Prepare XGBoost results for cross-classifier plotting.
stat_df <- formated_performance_report %>%
  rename(metrics=Metrics, mean=Value, group=ScanGroup) %>%
  mutate(model='XGBoost') %>% dplyr::select(-c("Dataset"))
# Run the sensitivity-analysis notebook first, then load the exported test-set metrics.
test_metrics <- read_csv(file.path(MODELING_TABLE_DIR, 'patient_metrics_split(smote_cv).csv'))
test_metrics <- test_metrics %>%
  dplyr::select(where(~ !any(is.na(.)))) %>%
  mutate(group = case_when(
    group == "Multi-center Manual" ~ "Manual",
    group == "H1/S1 Manual" ~ "DT_UIH_Manual",
    group == "H1/S1 DiscSeek" ~ "DT_UIH_AI",
    group == "H1/S1 DiscSeek-assisted" ~ "DT_UIH_Interactive",
    TRUE ~ group  # Keep other model group labels unchanged.
  ))
test_metrics <- test_metrics %>% pivot_longer(cols = -c(model, group), names_to = 'metrics', values_to = 'values')
test_metrics <- test_metrics %>%
  extract(
    values,
    into = c("mean", "ci_lower", "ci_upper"),
    regex = "(-?[0-9.]+)\\s*\\((-?[0-9.]+)\\s*[–-]\\s*(-?[0-9.]+)\\)",
    remove = FALSE,
    convert = TRUE
  )

test_metrics <- test_metrics %>% 
  bind_rows(stat_df) %>% 
  mutate(
    model = factor(model), 
    group = factor(group, 
                   levels = patient_scan_groups[c(4,3,1,2)], 
                   ordered = TRUE)
    )

p <- test_metrics %>% 
  filter(metrics %in% visible_metrics, model != 'nb') %>% 
  ggplot(aes(x=model, y=mean, fill = group)) + 
  geom_col(position = position_dodge(0.9)) + 
  geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), 
                position=position_dodge(0.9), 
                width=0.25) +
  facet_wrap(~ metrics, scales = 'free_y') + 
  labs(x='',y='')+
  theme_bw() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
    )

  print(p)
}
