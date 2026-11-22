# 02_descriptive.R
source("R/00_config.R")
library(ggplot2)
library(patchwork)
library(psych)
library(boot)
library(dplyr)
library(tidyr)
# Load the processed external validation dataset.
final_data <- read_rds(file.path(CLEAN_DATA_DIR, "processed_data.rds"))

# Convert categorical analysis variables to factors.
data <- final_data %>%
  mutate(across(c("ScanType", "ScanMachine", "ScanSite"), as.factor)) %>%
  mutate(across(contains("diagnosis"), as.factor)) %>%
  mutate(across(all_of(OUTLIER_CONFIG$sequence_params), as.factor)) %>%
  mutate(across(all_of(OUTLIER_CONFIG$discrete_params), as.factor))

# Assign the required levels to categorical and ordinal variables.
data <- data %>%  
  mutate(
    IVDType = factor(data$IVDType, levels = sort(unique(data$IVDType))),
    Sex = factor(Sex),
    Imaging.quality.score = factor(Imaging.quality.score, levels = sort(unique(data$Imaging.quality.score)))
  )
data$ScanMachine <- relevel(data$ScanMachine, ref="Philips")
data$ScanSite <- relevel(data$ScanSite, ref="QQHE")
data$ScanType<- relevel(data$ScanType, ref = "Manual")
data$Clinical.diagnosis <- relevel(factor(data$Clinical.diagnosis), ref = "ACS")
data$IVDType <- relevel(data$IVDType, ref = "C2C3")

# Construct overlapping planning groups for descriptive summaries.
data_cleaned <- data %>%
    mutate(ScanGroup = case_when(
    ScanType == "Manual" ~ "Multi-center Manual",
    ScanType == "AI" ~ "H1/S1 DiscSeek",
    ScanType == "Interactive" ~ "H1/S1 DiscSeek-assisted",
    TRUE ~ NA_character_
  ))

h1_manual <- data %>%
  filter(
    as.character(ScanSite) == "DT",
    as.character(ScanMachine) == "UIH",
    as.character(ScanType) == "Manual"
  ) %>%
  mutate(ScanGroup = "H1/S1 Manual")

data_cleaned <- bind_rows(data_cleaned, h1_manual) %>%
  filter(!is.na(ScanGroup))

# Supplementary Table S1 | Cohort size and age distribution ----
age_summary <- data_cleaned %>%
  filter(IVDType == "C2C3") %>%
  group_by(ScanGroup) %>%
  summarise(
    N = n(),
    Age_mean = mean(Age, na.rm = TRUE),
    Age_sd = sd(Age, na.rm = TRUE),
    Age_median = median(Age, na.rm = TRUE),
    Age_Q1 = quantile(Age, 0.25, na.rm = TRUE),
    Age_Q3 = quantile(Age, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Age_mean_sd = sprintf(
      "%.2f ± %.2f",
      Age_mean,
      Age_sd
    ),
    Age_median_IQR = sprintf(
      "%.0f (%.0f–%.0f)",
      Age_median,
      Age_Q1,
      Age_Q3
    )
  )

# Supplementary Table S1 | Categorical patient characteristics ----
# Summarise the within-group counts and proportions of a categorical variable.
summarise_category <- function(df, variable) {
  df %>%
    group_by(ScanGroup, .data[[variable]]) %>%
    summarise(
      n = n(),
      .groups = "drop"
    ) %>%
    group_by(ScanGroup) %>%
    mutate(
      total = sum(n),
      percentage = 100 * n / total,
      result = sprintf(
        "%d (%.2f%%)",
        n,
        percentage
      )
    ) %>%
    ungroup() %>%
    rename(Category = all_of(variable)) %>%
    mutate(Variable = variable)
}

# Gender
sex_summary <- summarise_category(
  data_cleaned %>% filter(IVDType=='C2C3'),
  "Sex"
)

# Cervical Spondylosis
diagnosis_summary <- summarise_category(
  data_cleaned %>% filter(IVDType=='C2C3'),
  "Clinical.diagnosis"
)

# Hospital-Vendor Setting
hospital_vendor_summary <- data_cleaned %>%
  filter(IVDType=='C2C3') %>% 
  mutate(
    Hospital = recode(
      as.character(ScanSite),
      "DT" = "H1",
      "QQHE" = "H2",
      "ZJ" = "H3"
    ),
    Vendor = recode(
      as.character(ScanMachine),
      "UIH" = "S1",
      "Siemens" = "S2",
      "Philips" = "S3",
      "GE" = "S4"
    )
  ) %>%
  count(
    ScanGroup,
    Hospital,
    Vendor,
    name = "n"
  ) %>%
  group_by(ScanGroup) %>%
  mutate(
    total = sum(n),
    percentage = 100 * n / total,
    result = sprintf(
      "%d (%.0f%%)",
      n,
      percentage
    )
  ) %>%
  ungroup()

# Supplementary Table S1 | Baseline characteristic comparisons ----
# Compare the complete Multi-center Manual cohort with both DiscSeek planning groups.
comparison1 <- data_cleaned %>% 
  filter(IVDType=="C2C3", ScanGroup != "H1/S1 Manual")%>%
  mutate(
    HospitalVendor = paste(
      ScanSite,
      ScanMachine,
      sep = "/"
    )
  )
# Compare the H1/S1 Manual subset with the non-H1/S1 Manual subset.
non_h1_manual <-  data_cleaned %>% 
  filter(ScanType == "Manual") %>% 
  anti_join(
    h1_manual %>% dplyr::select(ImageID),
    by = "ImageID")
comparison2 <- bind_rows(
  list(
    "H1/S1 Manual" = h1_manual,
    "Non-H1/S1 Manual" = non_h1_manual
  ),
  .id = "ScanGroup"
  ) %>% 
  filter(IVDType == "C2C3") %>%
  mutate(
    HospitalVendor = paste(
      ScanSite,
      ScanMachine,
      sep = "/"
    )
  )


# Compare age distributions.
age_test1 <- compare_means(
  Age ~ ScanGroup,
  data = comparison1,
  method = "kruskal.test"
)

age_test2 <- compare_means(
  Age ~ ScanGroup,
  data = comparison2,
  method = "kruskal.test"
)
# Compare categorical baseline characteristics.
categorical_test <- function(data, variable) {
  # Construct a contingency table.
  dis_table <- table(
    data[[variable]],
    data$ScanGroup
  )
  
  # Pearson Chi-square test
  chi_result <- suppressWarnings(
    chisq.test(dis_table)
  )
  
  # Return the test result.
  tibble(
    Parameter = variable,
    Chi_square_P = chi_result$p.value,
    Min_expected = min(chi_result$expected),
    Low_expected_percent =
      mean(chi_result$expected < 5) * 100
  )
}

categorical_results1 <- bind_rows(
  categorical_test(comparison1, "Sex"),
  categorical_test(comparison1, "Clinical.diagnosis")
)

categorical_results2 <- bind_rows(
  categorical_test(comparison2, "Sex"),
  categorical_test(comparison2, "Clinical.diagnosis")
)


# Figure 2k and Supplementary Table S3 | Image quality distributions ----
score_list <- summarise_category(
  data_cleaned %>% filter(IVDType=='C2C3'),
  "Imaging.quality.score"
) %>%
  rename(Imaging.quality.score=Category)

# Summarise the frequency of each observed image quality score by planning group.
score_list %>% 
  dplyr::select(c(ScanGroup, Imaging.quality.score, result, total)) %>%
  rename(IQ=result) %>%
  pivot_wider(names_from = Imaging.quality.score, values_from = c(IQ,total), names_sort = TRUE) %>%
  dplyr::select(-c(total_2, total_3, total_5)) %>%
  rename(Total=total_4)


# Table 2, Figure 2l, Tables S4–S5 | Radiological–clinical concordance ----
# Supplementary Table S5 | Bowker’s test of symmetry ----
diagnosis_levels <- c("ACS", "CSM", "CSR")
# Construct patient-level radiological–clinical diagnosis contingency tables.
dis_tables <- data_cleaned %>%
  distinct(ScanGroup, ImageID, .keep_all = TRUE) %>%
  group_by(ScanGroup) %>%
  group_modify(~ {
    dis_table <- table(
      Radiology = factor(
        .x$Radiology.diagnosis,
        levels = diagnosis_levels
      ),
      Clinical = factor(
        .x$Clinical.diagnosis,
        levels = diagnosis_levels
      )
    )
    
    as.data.frame(dis_table) %>%
      as_tibble()
  }) %>%
  ungroup()

# Apply Bowker’s test of symmetry within each planning group.
bowker_results <- dis_tables %>%
  group_by(ScanGroup) %>%
  group_modify(~ {
    dis_table <- xtabs(
      Freq ~ Radiology + Clinical,
      data = .x
    )
    
    test <- stats::mcnemar.test(dis_table)
    
    tibble(
      N = sum(dis_table),
      Statistic = unname(test$statistic),
      DF = unname(test$parameter),
      P.value = test$p.value
    )
  }) %>%
  ungroup()

# Combine the cell counts with Bowker test results.
result <- dis_tables %>%
  left_join(bowker_results, by = "ScanGroup")
result %>% pivot_wider(names_from = Clinical, values_from = c(Freq))

# Table 2 and Supplementary Table S4 | Cohen’s kappa ----
# Estimate kappa in each planning group and hospital–vendor strata for Figure 2l.
kappa_results <- data_cleaned %>%
  filter(IVDType == "C2C3") %>%
  filter(!(ScanGroup == "Multi-center Manual" & ScanMachine == "UIH")) %>%
  group_by(ScanGroup, ScanMachine, ScanSite) %>%
  summarise(
    kappa = cohen.kappa(cbind(Clinical.diagnosis, Radiology.diagnosis))$kappa,
    kappa_lower = cohen.kappa(cbind(Clinical.diagnosis, Radiology.diagnosis))$confid[1,1],
    kappa_upper = cohen.kappa(cbind(Clinical.diagnosis, Radiology.diagnosis))$confid[1,3],
    Scan.angle.offset = mean(Scan.angle.offset, na.rm = TRUE),
    Imaging.quality.score = mean(as.numeric(as.character(Imaging.quality.score)), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(ScanGroup = as.character(ScanGroup)) %>%
  dplyr::select(-c(ScanMachine, ScanSite)) %>%
  mutate(ScanGroup = ifelse(str_detect(ScanGroup, "H1/S1"), ScanGroup, "Manual(Single)"))
# Summarise agreement in the complete Multi-center Manual cohort.
kappa_results <- data_cleaned %>%
  filter(IVDType == "C2C3") %>% filter(str_detect(ScanGroup, "Multi-center")) %>% 
  summarise(kappa = cohen.kappa(cbind(Clinical.diagnosis, Radiology.diagnosis))$kappa,
            kappa_lower = cohen.kappa(cbind(Clinical.diagnosis, Radiology.diagnosis))$confid[1,1],
            kappa_upper = cohen.kappa(cbind(Clinical.diagnosis, Radiology.diagnosis))$confid[1,3],
            Scan.angle.offset = mean(Scan.angle.offset, na.rm = TRUE),
            Imaging.quality.score = mean(as.numeric(as.character(Imaging.quality.score)), na.rm = TRUE)) %>% 
  mutate(ScanGroup = "Multi-center Manual") %>% bind_rows(kappa_results)

# Supplementary Table S4 | Pairwise differences in kappa with bootstrap inference ----
calculate_kappa <- function(data, indices) {
  if (!missing(indices)) {
    data <- data[indices, ]
  }
  kappa_result <- cohen.kappa(data)
  return(kappa_result$kappa)
}

n_bootstraps <- 1000
patient_scan_groups <- c("H1/S1 DiscSeek", "H1/S1 DiscSeek-assisted", "H1/S1 Manual", "Multi-center Manual")
group_paires <- combn(patient_scan_groups, 2, simplify = FALSE)

pairwise_kappa_comp <- pmap_dfr(list(seq(1, length(group_paires))), function(i) {
  group1 <- group_paires[[i]][1]
  group2 <- group_paires[[i]][2]
  print(paste(group1, group2, sep = " VS "))
  if (all(str_detect(c(group1, group2), "Manual"))) {
    df_filter_1 <- data_cleaned %>% 
      filter(IVDType == "C2C3", ScanGroup == "Multi-center Manual", ScanMachine != "UIH") %>% 
      dplyr::select(c(Radiology.diagnosis, Clinical.diagnosis))
    df_filter_2 <- data_cleaned %>% 
      filter(IVDType == "C2C3", ScanGroup == "H1/S1 Manual") %>% 
      dplyr::select(c(Radiology.diagnosis, Clinical.diagnosis))
  } else {
    df_filter_1 <- data_cleaned %>% 
      filter(IVDType == "C2C3", str_detect(ScanGroup, group1)) %>% 
      dplyr::select(c(Radiology.diagnosis, Clinical.diagnosis))
    df_filter_2 <- data_cleaned %>% 
      filter(IVDType == "C2C3", str_detect(ScanGroup, group2)) %>% 
      dplyr::select(c(Radiology.diagnosis, Clinical.diagnosis))
  }
  
  kappa1_obs <- calculate_kappa(df_filter_1)
  kappa2_obs <- calculate_kappa(df_filter_2)
  delta_kappa_obs <- kappa1_obs - kappa2_obs
  
  boot_result1 <- boot(data = df_filter_1, statistic = calculate_kappa, R = n_bootstraps, strata = df_filter_1[[diag_var]])
  boot_result2 <- boot(data = df_filter_2, statistic = calculate_kappa, R = n_bootstraps, strata = df_filter_2[[diag_var]])
  delta_kappas <- boot_result1$t - boot_result2$t
  
  ci_result <- boot.ci(boot.out = list(t = delta_kappas, R = n_bootstraps, t0=delta_kappa_obs), 
                       type = "perc")  
  
  prop_greater <- mean(delta_kappas > 0, na.rm = TRUE)
  prop_less <- mean(delta_kappas < 0, na.rm = TRUE)
  p_value <- 2 * min(prop_greater, prop_less)
  tibble(
    Group1 = group1,
    Group2 = group2,
    n1 = nrow(df_filter_1),
    n2 = nrow(df_filter_2),
    kappa1_observed = kappa1_obs,
    kappa2_observed = kappa2_obs,
    delta_kappa_observed = delta_kappa_obs,
    delta_kappa_lower = ci_result$percent[4],
    delta_kappa_upper = ci_result$percent[5],
    p_value = p_value,
    n_bootstraps = n_bootstraps
  )
})

# Table 2 | Patient-level agreement rate and confidence interval ----
agreement_results <- data_cleaned %>%
  filter(IVDType == "C2C3") %>%
  distinct(ScanGroup, ImageID, .keep_all = TRUE) %>%
  filter(
    !is.na(Clinical.diagnosis),
    !is.na(Radiology.diagnosis)
  ) %>%
  group_by(ScanGroup) %>%
  group_modify(~ {
    
    # Construct a three-by-three diagnostic contingency table.
    dis_table <- table(
      factor(
        .x$Radiology.diagnosis,
        levels = diagnosis_levels
      ),
      factor(
        .x$Clinical.diagnosis,
        levels = diagnosis_levels
      )
    )
    
    # Count matching radiological and clinical diagnoses.
    x <- sum(diag(dis_table))
    
    # Count patient-level observations.
    n <- sum(dis_table)
    
    # Estimate agreement and its 95% Wilson score confidence interval.
    result <- stats::prop.test(
      x = x,
      n = n,
      correct = FALSE
    )
    
    tibble(
      N = n,
      Agreement_N = x,
      Agreement_rate = x / n,
      CI_lower = result$conf.int[1],
      CI_upper = result$conf.int[2]
    )
  }) %>%
  ungroup()

# Figure 2k–l | Image quality and diagnostic concordance plots ----
colors <- c("#3951A2", "#72AACF", "#CAE8F2", "#FEFBBA", "#FDB96B", "#EC5D3B", "#A80326")
color_palettes <- c(colors)[c(1,2,5,6)]
group_levels <- c("H1/S1 DiscSeek-assisted", "H1/S1 DiscSeek", "H1/S1 Manual", "Multi-center Manual") %>% rev()

# Figure 2k | Plot within-group image quality score proportions. ----
score_plot <- score_list %>% 
  dplyr::select(-c(result, Variable, total)) %>%
  complete(
    ScanGroup = group_levels, 
    Imaging.quality.score = levels(data_cleaned$Imaging.quality.score),
    fill = list(percentage = 0, n = 0)
  ) %>%
  mutate(ScanGroup = factor(ScanGroup, levels = group_levels),
         Imaging.quality.score = factor(Imaging.quality.score, levels = seq(5,2,by=-1)),
         ratio = percentage / 100) %>%
  mutate(
    ScanGroup = fct_relabel(ScanGroup, ~ {
      str_replace_all(., "H1/S1 ", "H1/S1\n") %>%
        str_replace_all(., "DiscSeek-", "DiscSeek-\n") %>%
        str_replace_all(., "Multi-center ", "Multi-center\n")
    })
  ) %>% 
  ggplot(aes(x = ScanGroup, y = ratio)) + 
  geom_bar(aes(group = Imaging.quality.score, fill = Imaging.quality.score), stat = "identity", width = 0.5, colour = "grey40", linewidth = 0.25) + 
  scale_fill_manual(values = setNames(color_palettes, seq(2,5))) + 
  labs(y = "Proportion", x = "Planning Group", fill = "IQ Score") + 
  guides(fill = guide_legend(reverse = TRUE)) + 
  scale_y_continuous(limits = c(0, 1.0), expand = expansion(mult = c(0, 0.25))) + 
  theme(
    # legend.title = element_text(size = 10, angle = 90, hjust = 0.5),
    legend.title.position = "left",
    legend.direction = "horizontal",
    legend.position = c(0.5, 0.9),
    legend.margin = ggplot2::margin(0,0,0,0),
    legend.spacing = unit(0, "pt"),
    legend.key.width = unit(12, "pt"),
    legend.key.height = unit(8, "pt"),
    # legend.background = element_rect(fill = "white", colour = "grey"),
    # legend.position = c(0.45, 1.06),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.box.spacing = unit(3,"pt"),
    legend.key = element_rect(linewidth = 0),
    axis.text = element_text(color = "black", size = 6),
    axis.text.x = element_text(size = 8),  #, angle = 12),
    axis.title = element_text(size = 10),
    panel.grid.minor = element_line(colour = "grey", linewidth=0.25),
    panel.grid.major = element_line(colour = "grey", linewidth = 0.25),
    panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
    # plot.margin = ggplot2::margin(t=25, unit="pt")
  )

# Figure 2l | Plot kappa estimates across planning groups and sites. ----
color_palettes <- colors[c(1,2,5,6,7)]
kappa_results <- kappa_results %>% 
  mutate(Group = factor(ScanGroup, levels = c("Multi-center Manual", "Manual(Single)", patient_scan_groups[c(3,1,2)]))) %>%
  mutate(Group = fct_relabel(
    Group, ~ {
      str_replace_all(
        ., c(
          "^H1/S1 DiscSeek-assisted$" = "H1/S1 DiscSeek-\nassisted"
        ))
    }))

kappa_list <- kappa_results %>% ggplot(aes(x = Imaging.quality.score, y = kappa)) + 
  geom_point(aes(fill = Group, size = Scan.angle.offset), colour = "grey40", stroke = 0.25, shape = 21, size = 5) + 
  scale_fill_manual(
    values = color_palettes, 
    labels = function(x) {
      x <- ifelse(grepl("Manual", x), "Manual", x)
      x <- str_remove(x, "H1/S1 ")
      x
      },
    name = "Planning\nMethod"
    ) + 
  guides(
    # color = "none", # guide_legend(order = 1),
    fill = guide_legend(order = 1, override.aes = list(size = 2)),
    size = guide_legend(order = 2, override.aes = list(stroke = 1))) + 
  # scale_size_continuous(
  #   breaks = c(2.7, 3.2, 3.7),
  #   range = c(2, 8),
  #   name = "Scan\nAngle\nOffset [°]") +
  labs(x = "IQ Score" , y = "Diagnostic Concordance") +
  stat_ellipse(
    data = kappa_results %>% 
      filter(Group %in% c("Manual(Single)","H1/S1 Manual")),
    type = "euclid",
    level = 0.1495,
    linewidth = 0.7,
    linetype = "dashed",
    colour = "#3951A2"
  ) +  
  geom_label(
    aes(label = Group), hjust = "right", nudge_x = -0.020, 
    data = kappa_results %>% 
      filter(Group %in% c("H1/S1 Manual", "H1/S1 DiscSeek-\nassisted", "H1/S1 DiscSeek")) %>% 
      mutate(
        Group = fct_relabel(Group, ~ {
          str_replace_all(., "DiscSeek-\nassisted", "DiscSeek-assisted")
        })
      ),
    show.legend = FALSE,
    size = 8,
    size.unit = "pt"
  ) + 
  geom_label(
    aes(label = Group), hjust = "left", nudge_x = 0.020, 
    data = kappa_results %>% filter(Group %in% c("Multi-center Manual")),
    show.legend = FALSE,
    size = 8,
    size.unit = "pt"
  ) + 
  scale_x_continuous(expand = expansion(mult = c(0.0,0.0)), limits = c(3.65, 4.15)) + 
  scale_y_continuous(expand = expansion(mult = c(0.0,0.0)), limits = c(0.45, 0.9)) + 
  # coord_fixed() +
  theme(
    # legend.title = element_text(size = 10, angle = 90, hjust = 0.5),
    legend.title.position = "top",
    # legend.direction = "horizontal",
    legend.position = "right",
    legend.margin = ggplot2::margin(0,0,0,0),
    legend.spacing = unit(5, "pt"),
    # legend.background = element_rect(fill = "white", colour = "grey"),
    # legend.position = c(0.45, 1.06),
    legend.title = element_text(size = 10, hjust = 0.5),
    legend.text = element_text(size = 8),
    legend.box.spacing = unit(3,"pt"),
    legend.key = element_rect(linewidth = 0),
    axis.text = element_text(color = "black", size = 6),
    axis.title = element_text(size = 10),
    panel.grid.minor = element_line(colour = "grey", linewidth=0.25),
    panel.grid.major = element_line(colour = "grey", linewidth = 0.25),
    panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
    # plot.margin = ggplot2::margin(t=25, unit="pt")
    )
ggsave(file.path(COMPARISON_FIGURE_DIR, "kappa_compare_1.png"), kappa_list, width = 3.5, height = 3.5, units = "in")  
  

# Combine Figure 2k and Figure 2l panels.
ggsave(file.path(COMPARISON_FIGURE_DIR, "kappa_score_compare.png"), 
       wrap_plots(list("score" = score_plot + 
                         labs(tag = "a") + 
                         theme(
                           plot.tag = element_text(
                             size = 10, 
                             margin = ggplot2::margin(0,0,0,0)
                             )
                           ),
                       "kappa" = kappa_list + 
                         labs(tag = "b", y = "Diagnostic Consistency") & 
                         theme(
                           plot.tag = element_text(
                             size = 10, 
                             margin = ggplot2::margin(0,0,0,0)
                             ),
                           legend.key.height = unit(18, "pt")
                           )
                       ),
                  widths = c(1.0,1.0), nrow = 1) + theme(plot.tag = element_text(size = 10, face = "bold")), 
       width = 180, height = 180*7.5/16, units = "mm", dpi = 300) 
