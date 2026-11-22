# DiscSeek | External scan angle offset analysis
# Outputs: Table 2, Supplementary Table S2, and Figure 2i–j.

# 1. Load and prepare external validation data ----
source("R/00_config.R")
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(patchwork)
library(lme4)
library(emmeans)
library(broom.mixed)
library(stringr)
library(purrr)
library(nortest)
# Load processed disc-level data.
data <- read_rds(file.path(CLEAN_DATA_DIR, "processed_data.rds"))

# Define factors, references, and mutually exclusive planning groups.
data <- data %>%
  mutate(across(c("ScanType", "ScanMachine", "ScanSite"), as.factor)) %>%
  mutate(across(contains("diagnosis"), as.factor)) %>%
  mutate(across(all_of(OUTLIER_CONFIG$sequence_params), as.factor)) %>%
  mutate(across(all_of(OUTLIER_CONFIG$discrete_params), as.factor))

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


data_cleaned <- data %>%
    unite("ScanGroup", ScanSite, ScanMachine, ScanType, remove = FALSE) %>%
  mutate(Age = (Age - min(data$Age)) / (max(data$Age) - min(data$Age)),
         ScanGroup = ifelse(str_detect(ScanGroup, "Manual") & !str_detect(ScanGroup, "UIH"), "Manual", ScanGroup)) %>%
  mutate(ScanGroup = str_replace(ScanGroup, "DT_UIH_", "H1/S1 "))

# Table 2 | Factors associated with axial scan angle offset ----
# Fit the multicenter manual mixed-effects model with disc level, hospital, vendor, and disc-level interactions.
model <- lmer(Scan.angle.offset ~ ScanSite*IVDType + ScanMachine*IVDType + (1|ImageID), 
              data = data_cleaned %>% filter(ScanType == "Manual"))
broom.mixed::tidy(car::Anova(model, type = "III"))
# Evaluate disc-level differences within the three H1/S1 platform-matched groups.
postho_df <- data_cleaned %>% 
  filter(ScanMachine == "UIH") %>%
  distinct(ScanGroup) %>%
  pmap_dfr(~ {
    current_group <- ..1
    df <- data_cleaned %>% filter(ScanGroup == current_group)
    compare_means(Scan.angle.offset~IVDType, df, method = "kruskal.test", p.adjust.method = "bonferroni") %>% 
      mutate(ScanGroup = current_group)
  })

# Table 2, Supplementary Table S2, Figure 2i–j | Scan angle offset comparisons ----
total_diff_posthoc_test_df = NULL
pairwise_comparisons <- NULL
adjusted_mean <- NULL
# H1/S1 platform-matched comparisons using distribution-based tests ----
# Assess the distribution of offsets within each planning group and cervical disc level.
min_p <-data_cleaned %>% 
  filter(ScanMachine == "UIH") %>%
  group_by(ScanGroup, IVDType) %>% 
  summarise(sample_size = n(), p = if (n() < 5000) {
    ifelse(n() > 3, shapiro.test(Scan.angle.offset)$p, NA)
  } else {
    ad.test(Scan.angle.offset)$p.value
  }, .groups = "drop") %>%
  ungroup()

# Assess offset distributions for each planning group across levels.
min_p <- data_cleaned %>% 
  filter(ScanMachine == "UIH") %>%
  group_by(ScanGroup) %>% 
  summarise(sample_size = n(), p = if (n() < 5000) {
    ifelse(n() > 3, shapiro.test(Scan.angle.offset)$p, NA)
  } else {
    ad.test(Scan.angle.offset)$p.value
  }, .groups = "drop") %>%
  ungroup() %>%
  mutate(IVDType = "Total") %>%
  bind_rows(min_p)
# Compare cervical disc levels within H1/S1 planning groups.
postho_df <- data_cleaned %>% 
  filter(ScanMachine == "UIH") %>%
  distinct(ScanGroup) %>%
  pmap_dfr(~ {
    current_group <- ..1
    df <- data_cleaned %>% filter(ScanGroup == current_group)
    compare_means(Scan.angle.offset~IVDType, df, method = "wilcox.test", p.adjust.method = "bonferroni") %>% 
      mutate(ScanGroup = current_group)
  })
one_way_test_df <- postho_df %>%
  distinct(ScanGroup, group1, group2) %>%
  pmap_dfr(~ {
    current_group <- ..1
    ivd_type1 <- ..2
    ivd_type2 <- ..3
    p <- min_p %>% filter((IVDType==ivd_type1)|(IVDType == ivd_type2), ScanGroup == current_group) %>% pull(p) %>% min()
    print(paste(c(current_group, ivd_type1, ivd_type2, p)))
    data_filtered <- data_cleaned %>%
      filter(ScanGroup == current_group, (IVDType == ivd_type1)|(IVDType == ivd_type2))
    if (length(unique(data_filtered$IVDType)) > 2) {
      # Handle comparisons involving more than two eligible levels.
      pairwise = TRUE
    } else{
      pairwise = FALSE
    }
    test_result <- perform_diff_test(data_filtered,
                                     as.formula("Scan.angle.offset ~ IVDType"),
                                     group_var = "IVDType",
                                     is_norm = p > 0.05,
                                     paired = FALSE,
                                     pairwise = pairwise) %>%
      mutate(ScanGroup = current_group, group1 = ivd_type1, group2 = ivd_type2)
    test_result
  })
total_diff_posthoc_test_df <- bind_rows(
  total_diff_posthoc_test_df, 
  one_way_test_df %>% 
    unite("IVDType", c("group1","group2"), sep = " - ")
  )

# Compare independent planning groups at each cervical disc level.
postho_df <- compare_means(Scan.angle.offset~ScanGroup, data_cleaned %>% filter(ScanMachine == "UIH"), 
                           group.by = "IVDType", method = "wilcox.test", p.adjust.method = "bonferroni") %>%
  mutate(IVDType = as.character(IVDType))
one_way_test_df <- postho_df %>%
  distinct(IVDType, group1, group2) %>%
  pmap_dfr(~ {
    current_group <- ..1
    type1 <- ..2
    type2 <- ..3
    p <- min_p %>% filter((ScanGroup==type1)|(ScanGroup == type2), IVDType == current_group) %>% pull(p) %>% min()
    print(paste(c(current_group, type1, type2, p)))
    data_filtered <- data_cleaned %>%
      filter(IVDType == current_group, (ScanGroup == type1)|(ScanGroup == type2))
    if (length(unique(data_filtered$ScanGroup)) > 2) {
      # Handle pairwise comparisons between two planning groups.
      pairwise = TRUE
    } else{
      pairwise = FALSE
    }
    test_result <- perform_diff_test(data_filtered,
                                     as.formula("Scan.angle.offset ~ ScanGroup"),
                                     group_var = "ScanGroup",
                                     is_norm = p > 0.05,
                                     paired = FALSE,
                                     pairwise = pairwise) %>%
      mutate(ScanGroup = current_group, group1 = type1, group2 = type2)
    test_result
  })
total_diff_posthoc_test_df <- bind_rows(total_diff_posthoc_test_df, 
                                        one_way_test_df %>% 
                                          mutate(ScanGroup = paste0("(",group1,") - (",group2, ")")) %>% 
                                          dplyr::select(-c(group1, group2)))
# Compare independent planning groups across disc levels.
one_way_test_df <- postho_df %>%
  distinct(group1, group2) %>%
  pmap_dfr(~ {
    type1 <- ..1
    type2 <- ..2
    p <- min_p %>% filter((ScanGroup==type1)|(ScanGroup == type2), IVDType == "Total") %>% pull(p) %>% min()
    print(paste(c(type1, type2, p)))
    data_filtered <- data_cleaned %>%
      filter((ScanGroup == type1)|(ScanGroup == type2))
    if (length(unique(data_filtered$ScanGroup))> 2) {
      # Handle comparisons of two planning groups.
      pairwise = TRUE
    } else{
      pairwise = FALSE
    }
    test_result <- perform_diff_test(data_filtered,
                                     as.formula("Scan.angle.offset ~ ScanGroup"),
                                     group_var = "ScanGroup",
                                     is_norm = p > 0.05,
                                     paired = FALSE,
                                     pairwise = pairwise) %>%
      mutate(IVDType = "Total", group1 = type1, group2 = type2)
    test_result
  })
total_diff_posthoc_test_df <- bind_rows(total_diff_posthoc_test_df, 
                                        one_way_test_df %>% 
                                          mutate(ScanGroup = paste0("(",group1,") - (",group2, ")")) %>% 
                                          dplyr::select(-c(group1, group2)))

# Multicenter manual comparisons using mixed-effects models and estimated marginal means ----
# Estimate disc-level marginal means and within-group pairwise contrasts for multicenter manual planning.
model <- lmer(Scan.angle.offset ~ IVDType + (1|ScanSite) + (1|ScanMachine) + (1|ImageID), data = data_cleaned %>% filter(ScanType == "Manual"))
emmeans_results_1 <- emmeans(model, pairwise ~ IVDType, adjust = "bonferroni", pbkrtest.limit = 6000)
gc()
pairwise_comparisons_1 <- summary(emmeans_results_1$contrasts)
pairwise_comparisons <- bind_rows(pairwise_comparisons, 
                                  pairwise_comparisons_1 %>% 
                                    mutate(ScanGroup = "Manual") %>%
                                    rename(IVDType = contrast)
                                  )
adjusted_mean_1 <- emmeans_results_1$emmeans %>%
  as.data.frame() %>% 
  mutate(ScanGroup = "Manual")
adjusted_mean <- bind_rows(adjusted_mean, adjusted_mean_1)

# Compare multicenter manual planning with H1/S1 groups at individual cervical disc levels.
# Compare all manual examinations with H1/S1 DiscSeek planning; use ScanType to distinguish approaches.
model <- lmer(Scan.angle.offset ~ IVDType*ScanType + (1|ScanSite) + (1|ScanMachine) + (1|ImageID), data = data_cleaned %>% filter(ScanType != "Interactive"))
emmeans_results_2 <- emmeans(model, pairwise ~ ScanType | IVDType, adjust = "bonferroni", pbkrtest.limit = 9000)
gc()
pairwise_comparisons_2 <- summary(emmeans_results_2$contrasts) %>% 
  rename(ScanGroup = contrast)
# Compare all manual examinations with H1/S1 DiscSeek-assisted planning; use ScanType to distinguish approaches.
model <- lmer(Scan.angle.offset ~ IVDType*ScanType + (1|ScanSite) + (1|ScanMachine) + (1|ImageID), data = data_cleaned %>% filter(ScanType != "AI"))
emmeans_results_3 <- emmeans(model, pairwise ~ ScanType | IVDType, adjust = "bonferroni", pbkrtest.limit = 9000)
gc()
pairwise_comparisons_2 <- summary(emmeans_results_3$contrasts) %>% 
  rename(ScanGroup = contrast) %>% bind_rows(pairwise_comparisons_2)
# Compare H1/S1 Manual with the independent non-H1/S1 Manual subset; use ScanGroup.
model <- lmer(Scan.angle.offset ~ IVDType*ScanGroup + (1|ScanSite) + (1|ScanMachine) + (1|ImageID), data = data_cleaned %>% filter(ScanType == "Manual"))
emmeans_results_4 <- emmeans(model, pairwise ~ ScanGroup | IVDType, adjust = "bonferroni", pbkrtest.limit = 9000)
gc()
pairwise_comparisons_2 <- summary(emmeans_results_4$contrasts) %>% 
  rename(ScanGroup = contrast) %>% bind_rows(pairwise_comparisons_2) %>%
  mutate(ScanGroup = str_replace_all(
    ScanGroup, c(
      "AI" = "H1/S1 Auto",
      "Interactive" = "H1/S1 Interactive"
    )))
pairwise_comparisons <- bind_rows(pairwise_comparisons, pairwise_comparisons_2)

# Estimate marginal planning-group contrasts averaged across disc levels.
model <- lmer(Scan.angle.offset ~ ScanType + (1|ScanSite) + (1|ScanMachine) + (1|ImageID), data = data_cleaned %>% filter(ScanType != "Interactive"))
emmeans_results_5 <- emmeans(model, pairwise ~ ScanType, adjust = "bonferroni", pbkrtest.limit = 9000)
gc()
pairwise_comparisons_3 <- summary(emmeans_results_5$contrasts) %>% 
  rename(ScanGroup = contrast) %>%
  mutate(ScanGroup = str_replace_all(
    ScanGroup, c(
      "AI" = "H1/S1 AI",
      "Interactive" = "H1/S1 Interactive"
    )),
    IVDType = "Total")
model <- lmer(Scan.angle.offset ~ ScanType + (1|ScanSite) + (1|ScanMachine) + (1|ImageID), data = data_cleaned %>% filter(ScanType != "AI"))
emmeans_results_6 <- emmeans(model, pairwise ~ ScanType, adjust = "bonferroni", pbkrtest.limit = 9000)
gc()
pairwise_comparisons_3 <- summary(emmeans_results_6$contrasts) %>% 
  rename(ScanGroup = contrast) %>% 
  mutate(IVDType='Total') %>%
  mutate(ScanGroup = str_replace_all(
    ScanGroup, c(
      "AI" = "H1/S1 Auto",
      "Interactive" = "H1/S1 Interactive"
    ))) %>%
  bind_rows(pairwise_comparisons_3)

model <- lmer(Scan.angle.offset ~ ScanGroup + (1|ScanSite) + (1|ScanMachine) + (1|ImageID), data = data_cleaned %>% filter(ScanType == "Manual"))
emmeans_results_7 <- emmeans(model, pairwise ~ ScanGroup, adjust = "bonferroni", pbkrtest.limit = 9000)
gc()
pairwise_comparisons_3 <- summary(emmeans_results_6$contrasts) %>% 
  rename(ScanGroup = contrast) %>% 
  mutate(IVDType='Total') %>%
  bind_rows(pairwise_comparisons_3)
pairwise_comparisons <- bind_rows(pairwise_comparisons, pairwise_comparisons_3)

# Optional analysis | Alternative mixed-effects model for H1/S1 planning groups ----
# Optional: this platform-matched model is not required for the principal comparisons.
# model <- lmer(Scan.angle.offset ~ IVDType * ScanGroup + (1 | ImageID),
#               data = data_cleaned %>% filter(ScanMachine == "UIH"))
# emmeans_results_8 <- emmeans(model, pairwise ~ IVDType | ScanGroup,
#                              adjust = "bonferroni", pbkrtest.limit = 6000)
# adjusted_mean_2 <- as.data.frame(emmeans_results_8$emmeans)
# adjusted_mean <- bind_rows(adjusted_mean, adjusted_mean_2)
# pairwise_comparisons_4 <- summary(emmeans_results_8$contrasts) %>%
#   rename(IVDType = contrast)
# pairwise_comparisons <- bind_rows(pairwise_comparisons, pairwise_comparisons_4)
# emmeans_results_9 <- emmeans(model, pairwise ~ ScanGroup | IVDType,
#                              adjust = "bonferroni", pbkrtest.limit = 6000)
# pairwise_comparisons_5 <- summary(emmeans_results_9$contrasts) %>%
#   rename(ScanGroup = contrast)
# pairwise_comparisons <- bind_rows(pairwise_comparisons, pairwise_comparisons_5)
# emmeans_results_10 <- emmeans(model, pairwise ~ ScanGroup,
#                               adjust = "bonferroni", pbkrtest.limit = 6000)
# adjusted_mean_3 <- as.data.frame(emmeans_results_10$emmeans) %>%
#   mutate(IVDType = "Total")
# mean_results <- bind_rows(adjusted_mean, adjusted_mean_3)

# Combine observed means, standard deviations, and model-adjusted estimates.
mean_results <- bind_rows(data_cleaned %>% 
                            filter(ScanMachine == "UIH") %>%
                            group_by(ScanGroup, IVDType) %>% 
                            summarise(mean = mean(Scan.angle.offset), sd = sd(Scan.angle.offset)),
                          data_cleaned %>% 
                            filter(ScanMachine == "UIH") %>%
                            group_by(ScanGroup) %>% 
                            summarise(mean = mean(Scan.angle.offset), sd = sd(Scan.angle.offset)) %>%
                            mutate(IVDType = "Total")
                          )
mean_results <- data_cleaned %>% 
  filter(ScanType == "Manual") %>% 
  group_by(ScanType, IVDType) %>% 
  summarise(mean=mean(Scan.angle.offset), 
            sd=sd(Scan.angle.offset)) %>% 
  rename(ScanGroup=ScanType) %>% 
  bind_rows(mean_results)

mean_results <- data_cleaned %>% filter(ScanType == "Manual") %>% 
  group_by(ScanType) %>% 
  summarise(mean=mean(Scan.angle.offset), 
            sd=sd(Scan.angle.offset)) %>% 
  rename(ScanGroup=ScanType) %>%
  mutate(IVDType = "Total") %>%
  bind_rows(mean_results)

mean_results <- mean_results%>%
  left_join(adjusted_mean, by=c("ScanGroup", "IVDType"))

# Export estimated marginal means and pairwise comparisons.
write_csv(mean_results, file.path(COMPARISON_TABLE_DIR, "angle_offset_stat.csv"))
write_csv(pairwise_comparisons, file.path(COMPARISON_TABLE_DIR, "angle_offset_diff_posthoc_test_ind_lm.csv"))


# Figure 2i–j | Prepare scan angle offset visualizations ----
# Apply manuscript labels to planning groups and cervical disc levels.
ivd_order_levels <- c("C2C3", "C3C4", "C4C5", "C5C6", "C6C7", "C7T1", "Total")
ivd_level_names <- c("C2/C3", "C3/C4", "C4/C5", "C5/C6", "C6/C7", "C7/T1", "Total")
ivd_level_dict <- setNames(ivd_level_names, ivd_order_levels)
group_level_names <- c("Multi-center\nManual", "H1/S1\nManual", "H1/S1\nDiscSeek", "H1/S1\nDiscSeek-\nassisted")
# Configure plot colors.
color_palettes <- c("#3951A2", "#72AACF", "#CAE8F2", "#FEFBBA", "#FDB96B", "#EC5D3B", "#A80326")
# Read saved marginal estimates and standardize group labels for plotting.
adjusted_mean <- read_csv(file.path(COMPARISON_TABLE_DIR, "angle_offset_stat.csv"))
pairwise_comparisons <- read_csv(file.path(COMPARISON_TABLE_DIR, "angle_offset_diff_posthoc_test_ind_lm.csv"))
adjusted_mean <- adjusted_mean %>% mutate(ScanGroup = str_replace_all(ScanGroup, c(
  "^H1/S1 AI$" = "H1/S1\nDiscSeek",          
  "^H1/S1 Interactive$" = "H1/S1\nDiscSeek-\nassisted", 
  "^H1/S1 Manual$" = "H1/S1\nManual",      
  "^Manual$" = "Multi-center\nManual"
  )),
  IVDType = factor(IVDType, levels = ivd_order_levels)) %>% 
  mutate(ScanGroup = factor(ScanGroup, levels = group_level_names),
         IVDType = fct_relabel(IVDType, ~{
           str_replace_all(., ivd_level_dict)}))
pairwise_comparisons <- pairwise_comparisons %>% 
  mutate(p.signif = case_when(p.value < 0.0001 ~ "****", 
                              p.value < 0.001 ~ "***", 
                              p.value < 0.01 ~ "**", 
                              p.value < 0.05 ~ "*", 
                              TRUE ~ "ns"))

# Combine nonparametric tests and mixed-effects model post hoc results.
pairwise_comparisons <- anti_join(pairwise_comparisons, 
          total_diff_posthoc_test_df, 
          by = c("IVDType", "ScanGroup")) %>% 
  mutate(method = "linear mixed effect model", 
         p.format = p_format(p.value, digits = 2)) %>% 
  bind_rows(total_diff_posthoc_test_df %>% 
              rename(p.value = p) %>% 
              mutate(p.format = p_format(p.value, digits = 2)) %>% 
              dplyr::select(-c(.y.,p.adj)))


# Figure 2j | Within-group variation across cervical disc levels ----
stat_df <- adjusted_mean %>% filter(IVDType != "Total")
# Draw the disc-level comparison plot.
p1 <- stat_df %>% ggplot(aes(x = ScanGroup, y = emmean, group = IVDType)) + 
  geom_bar(aes(fill = IVDType), stat = "identity", position = position_dodge2(0.7, preserve = "single"), 
         width = 0.7, colour = "grey40", linewidth = 0.25) +
  labs(x = "Planning Group", y = "Axial Scan Angle Offset [°]", fill = "Disc Level") + 
  scale_fill_manual(values = color_palettes) + 
  # guides(fill = guide_legend(nrow = 1)) + 
  theme(legend.position = "top", legend.margin = ggplot2::margin(b=0),
        legend.title = element_text(size = 10), legend.key.height = unit(5, "pt"),
        legend.spacing = unit(3,"pt"), legend.text = element_text(size = 8), legend.box.spacing = unit(0,"mm"),
        legend.key = element_rect(linewidth = 0, color = "white"), axis.text = element_text(color = "black", size = 8),
        axis.title = element_text(size = 10), axis.text.y = element_text(size = 6),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.25),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.25),
        panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 8, margin = ggplot2::margin(t=0,b=1)))

gg_data <- ggplot_build(p1)$data[[1]]
x_pos <- list()
for (scan_group in levels(adjusted_mean$ScanGroup)) {
  x_pos[[scan_group]] <- gg_data %>% 
    filter((x < match(scan_group, levels(adjusted_mean$ScanGroup)) + 0.5) & (x > match(scan_group, levels(adjusted_mean$ScanGroup)) - 0.5)) %>% 
    pull(x)%>%unique() %>% sort()
}

# Extract within-group post hoc test results.
posthoc_df <- pairwise_comparisons %>% filter(str_detect(IVDType, "-")) %>%
  separate("IVDType", c("group1", "group2"), sep = "-") %>%
  mutate(group1 = str_remove_all(group1, "\\(|\\)| "),
         group2 = str_remove_all(group2, "\\(|\\)| ")) %>% 
  mutate(group1 = str_remove_all(group1, " "),
         group2 = str_remove_all(group2, " ")) %>% 
  mutate(ScanGroup = str_replace_all(ScanGroup, c(
  "^H1/S1 AI$" = "H1/S1\nDiscSeek",          
  "^H1/S1 Interactive$" = "H1/S1\nDiscSeek-\nassisted", 
  "^H1/S1 Manual$" = "H1/S1\nManual",      
  "^Manual$" = "Multi-center\nManual"
  )),
  group1 = factor(group1, levels = ivd_order_levels),
  group2 = factor(group2, levels = ivd_order_levels)) %>% 
  mutate(ScanGroup = factor(ScanGroup, levels = c("Multi-center\nManual", "H1/S1\nManual", "H1/S1\nDiscSeek", "H1/S1\nDiscSeek-\nassisted")),
         group1 = fct_relabel(group1, ~{
           str_replace_all(., ivd_level_dict)}),
         group2 = fct_relabel(group2, ~{
           str_replace_all(., ivd_level_dict)}))
max_y <- stat_df %>% group_by(ScanGroup) %>% summarise(max_y = max(emmean + sd), .groups = "drop")
y_step <- max(max_y$max_y)

# Calculate annotation positions for statistical significance.
posthoc_df <- posthoc_df %>% 
  left_join(max_y, by = c("ScanGroup")) %>%
  mutate(ScanGroup = factor(ScanGroup, levels = levels(stat_df$ScanGroup)))
posthoc_df <- posthoc_df %>%
  group_by(ScanGroup) %>% 
  mutate(x = map(ScanGroup, ~ x_pos[[.]]),
         group = list(sort(unique(c(group1, group2)), 
                           method = "radix")[order(match(sort(unique(c(group1, group2)), 
                                                              method = "radix"), 
                                                         group_level_names))]),
         id1 = match(group1, group[[1]]),
         id2 = match(group2, group[[1]]),
         x_1 = map2_dbl(x, id1, ~ .x[.y]),
         x_2 = map2_dbl(x, id2, ~ .x[.y])) %>% ungroup()
posthoc_df <- posthoc_df %>% 
  filter(p.signif != "ns") %>%
  group_by(ScanGroup) %>%
  mutate(
    x_min = pmin(x_1, x_2),
    x_max = pmax(x_1, x_2),
    x_text = (x_min + x_max) / 2,
    y = max_y + (row_number()) * 0.05 * y_step, 
  ) %>% ungroup() %>%
  dplyr::select(-c(x, id1, id2, x_1, x_2, max_y, group)) %>%
  arrange(ScanGroup, group1, group2)
posthoc_df[9,"y"] = posthoc_df[2, "y"]
posthoc_df[8,"y"] = posthoc_df[1, "y"]

p1 <- p1 +
  geom_segment(aes(x = x_min, xend = x_max, y = y, yend = y), 
               data = posthoc_df , inherit.aes = FALSE, linewidth = 0.25) + 
  geom_segment(aes(x = x_min, xend = x_min, y = 0.99 * y, yend = y), 
               data = posthoc_df , inherit.aes = FALSE, linewidth = 0.25) + 
  geom_segment(aes(x = x_max, xend = x_max, y = 0.99 * y, yend = y), 
               data = posthoc_df , inherit.aes = FALSE, linewidth = 0.25) + 
  geom_text(aes(label = p.signif, x = x_text, y = y),
            data = posthoc_df, inherit.aes = FALSE, size = 3, vjust = 0, nudge_y = -0.05)

# Calculate standard deviation whisker positions.
error_df <- stat_df %>% group_by(ScanGroup) %>%
  mutate(x.pos = map(ScanGroup, ~ x_pos[[.]]),  
         group = list(sort(unique(IVDType), 
                           method = "radix")[order(match(sort(unique(IVDType), 
                                                              method = "radix"), 
                                                         group_level_names))]),
         id = match(IVDType, group[[1]]),
         x = map2_dbl(x.pos, id, ~ .x[.y])) %>% 
  ungroup() %>% 
  dplyr::select(-c(x.pos, group, id))
bar_width <- min(gg_data$xmax - gg_data$xmin)

p1 <- p1 +
  geom_segment(aes(x = x, xend = x, y = mean - sd, yend = mean + sd), 
               data = error_df , inherit.aes = FALSE, linewidth = 0.25) + 
  geom_segment(aes(x = x - bar_width * 0.1, xend = x + bar_width * 0.1, y = mean - sd, yend = mean - sd), 
               data = error_df , inherit.aes = FALSE, linewidth = 0.25) + 
  geom_segment(aes(x = x - bar_width * 0.1, xend = x + bar_width * 0.1, y = mean + sd, yend = mean + sd), 
               data = error_df , inherit.aes = FALSE, linewidth = 0.25) 

p1 <- p1 + scale_y_continuous(
  limits = c(0, NA),
  expand = expansion(mult = c(0, 0.1))
)

# Optional visualization | Overall planning-group differences already represented in Figure 2i ----
# max_y <- max(posthoc_df$y)
# posthoc_df <- pairwise_comparisons %>% filter(str_detect(ScanGroup, "-"), IVDType == "Total") %>%
#   separate("ScanGroup", c("group1", "group2"), sep = "-") %>%
#   mutate(group1 = str_remove_all(str_trim(group1), "\\(|\\)"),
#          group2 = str_remove_all(str_trim(group2), "\\(|\\)")) %>%
#   mutate(group1 = str_replace_all(group1, c(
#     "H1/S1 AI" = "H1/S1\nDiscSeek",          
#     "H1/S1 Interactive" = "H1/S1\nDiscSeek-\nassisted", 
#     "H1/S1 Manual" = "H1/S1\nManual",      
#     "^Manual$" = "Multi-center\nManual"
#   )),
#   group2 = str_replace_all(group2, c(
#     "H1/S1 AI" = "H1/S1\nDiscSeek",          
#     "H1/S1 Interactive" = "H1/S1\nDiscSeek-\nassisted", 
#     "H1/S1 Manual" = "H1/S1\nManual",      
#     "^Manual$" = "Multi-center\nManual"
#   ))) %>% 
#   mutate(group1 = factor(group1, levels = group_level_names),
#          group2 = factor(group2, levels = group_level_names),
#          max_y = max_y)
# 
# posthoc_df <- posthoc_df %>%
#   mutate(x_1 = map_dbl(group1, ~ mean(x_pos[[.x]])),
#          x_2 = map_dbl(group2, ~ mean(x_pos[[.x]])))
# 
# posthoc_df <- posthoc_df %>% 
#   filter(p.signif != "ns") %>%
#   arrange(desc(group1), desc(group2)) %>% 
#   mutate(
#     x_min = pmin(x_1, x_2),
#     x_max = pmax(x_1, x_2),
#     x_text = (x_min + x_max) / 2,
#     y = max_y + (row_number()-2) * 0.05 * y_step, 
#   ) %>% 
#   dplyr::select(-c(x_1, x_2, max_y))
# 
# p1_sig <- p1 +
#   coord_cartesian(ylim = NULL)
# p1_sig <- p1_sig + 
#   geom_segment(aes(x = x_min, xend = x_max, y = y, yend = y), 
#                data = posthoc_df , inherit.aes = FALSE, linewidth = 0.25) + 
#   geom_segment(aes(x = x_min, xend = x_min, y = 0.98 * y, yend = y), 
#                data = posthoc_df , inherit.aes = FALSE, linewidth = 0.25) + 
#   geom_segment(aes(x = x_max, xend = x_max, y = 0.98 * y, yend = y), 
#                data = posthoc_df , inherit.aes = FALSE, linewidth = 0.25) + 
#   geom_text(aes(label = p.signif, x = x_text, y = y),
#             data = posthoc_df, inherit.aes = FALSE, size = 3, vjust = 0, nudge_y = -0.05)
# 
# p1_sig <- p1_sig + scale_y_continuous(
#   limits = c(0, NA),
#   expand = expansion(mult = c(0, 0.1))
# )
ggsave(file.path(COMPARISON_FIGURE_DIR, "intra_angle_diff_test_by_ivdtype.png"), p1_sig, width = 180, height = 180*7/16, units = "mm", dpi = 300)

# Figure 2i | Between-group offsets at each cervical disc level ----
stat_df <- adjusted_mean
p2 <- stat_df %>% ggplot(aes(group = ScanGroup, y = emmean, x = IVDType)) + 
  geom_bar(aes(fill = ScanGroup), stat = "identity", position = position_dodge2(0.7, preserve = "single"), 
           width = 0.7, colour = "grey40", linewidth = 0.25) +
  labs(fill = "Planning Group", y = "Axial Scan Angle Offset [°]", x = "Disc Level") + 
  scale_fill_manual(values = color_palettes[c(1,2,5,6)]) + 
  guides(fill = guide_legend(nrow = 1)) + 
  theme(legend.position = "top", legend.margin = ggplot2::margin(b=0),
        legend.title = element_text(size = 10), legend.key.height = unit(5, "pt"),
        legend.spacing = unit(3,"pt"), legend.text = element_text(size = 8), legend.box.spacing = unit(0,"mm"),
        legend.key = element_rect(linewidth = 0, color = "white"), axis.text = element_text(color = "black", size = 8),
        axis.title = element_text(size = 10), axis.text.y = element_text(size = 6),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.25),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.25),
        panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 8, margin = ggplot2::margin(t=0,b=1)))

gg_data <- ggplot_build(p2)$data[[1]]
x_pos <- list()
for (ivd_type in levels(stat_df$IVDType)) {
  x_pos[[ivd_type]] <- gg_data %>% 
    filter((x < match(ivd_type, ivd_level_names) + 0.5) & (x > match(ivd_type, ivd_level_names) - 0.5)) %>% 
    pull(x)%>%unique() %>% sort()
}

# Prepare between-group pairwise significance annotations.
posthoc_df <- pairwise_comparisons %>% filter(str_detect(ScanGroup, "-")) %>%
  separate("ScanGroup", c("group1", "group2"), sep = "-") %>%
  mutate(group1 = str_remove_all(str_trim(group1), "\\(|\\)"),
         group2 = str_remove_all(str_trim(group2), "\\(|\\)")) %>%
  mutate(group1 = str_replace_all(group1, c(
    "H1/S1 AI" = "H1/S1\nDiscSeek",          
    "H1/S1 Interactive" = "H1/S1\nDiscSeek-\nassisted", 
    "H1/S1 Manual" = "H1/S1\nManual",      
    "^Manual$" = "Multi-center\nManual"
  )),
  group2 = str_replace_all(group2, c(
    "H1/S1 AI" = "H1/S1\nDiscSeek",          
    "H1/S1 Interactive" = "H1/S1\nDiscSeek-\nassisted", 
    "H1/S1 Manual" = "H1/S1\nManual",      
    "^Manual$" = "Multi-center\nManual"
  ))) %>% 
  mutate(group1 = factor(group1, levels = group_level_names),
         group2 = factor(group2, levels = group_level_names),
         IVDType = factor(IVDType, levels = ivd_order_levels)) %>% 
  mutate(IVDType = fct_relabel(IVDType, ~{
           str_replace_all(., ivd_level_dict)}))
max_y <- stat_df %>% group_by(IVDType) %>% summarise(max_y = max(emmean + sd), .groups = "drop")
y_step <- max(max_y$max_y)

posthoc_df <- posthoc_df %>% 
  left_join(max_y, by = c("IVDType")) %>%
  mutate(IVDType = factor(IVDType, levels = ivd_level_names))
posthoc_df <- posthoc_df %>%
  group_by(IVDType) %>% 
  mutate(x = map(IVDType, ~ x_pos[[.]]),
         group = list(sort(unique(c(group1, group2)), 
                           method = "radix")[order(match(sort(unique(c(group1, group2)), 
                                                              method = "radix"), 
                                                         ivd_level_names))]),
         id1 = match(group1, group[[1]]),
         id2 = match(group2, group[[1]]),
         x_1 = map2_dbl(x, id1, ~ .x[.y]),
         x_2 = map2_dbl(x, id2, ~ .x[.y])) %>% ungroup()
posthoc_df <- posthoc_df %>% 
  filter(p.signif != "ns") %>%
  group_by(IVDType) %>%
  mutate(
    x_min = pmin(x_1, x_2),
    x_max = pmax(x_1, x_2),
    x_text = (x_min + x_max) / 2,
    y = max_y + (row_number()) * 0.1 * y_step, 
  ) %>% ungroup() %>%
  dplyr::select(-c(x, id1, id2, x_1, x_2, max_y, group)) %>%
  arrange(IVDType, group1, group2)

p2 <- p2 +
  geom_segment(aes(x = x_min, xend = x_max, y = y, yend = y), 
               data = posthoc_df , inherit.aes = FALSE, linewidth = 0.25) + 
  geom_segment(aes(x = x_min, xend = x_min, y = 0.99 * y, yend = y), 
               data = posthoc_df , inherit.aes = FALSE, linewidth = 0.25) + 
  geom_segment(aes(x = x_max, xend = x_max, y = 0.99 * y, yend = y), 
               data = posthoc_df , inherit.aes = FALSE, linewidth = 0.25) + 
  geom_text(aes(label = p.signif, x = x_text, y = y),
            data = posthoc_df, inherit.aes = FALSE, size = 3, vjust = 0, nudge_y = -0.2)

error_df <- stat_df %>% group_by(IVDType) %>%
  mutate(x.pos = map(IVDType, ~ x_pos[[.]]),
         group = list(sort(unique(ScanGroup), 
                           method = "radix")[order(match(sort(unique(ScanGroup), 
                                                              method = "radix"), 
                                                         ivd_level_names))]),
         id = match(ScanGroup, group[[1]]),
         x = map2_dbl(x.pos, id, ~ .x[.y])) %>% 
  ungroup() %>% 
  dplyr::select(-c(x.pos, group, id)) %>%
  arrange(IVDType, ScanGroup)
bar_width <- min(gg_data$xmax - gg_data$xmin)

p2 <- p2 +
  geom_segment(aes(x = x, xend = x, y = mean - sd, yend = mean + sd), 
               data = error_df , inherit.aes = FALSE, linewidth = 0.25) + 
  geom_segment(aes(x = x - bar_width * 0.1, xend = x + bar_width * 0.1, y = mean - sd, yend = mean - sd), 
               data = error_df , inherit.aes = FALSE, linewidth = 0.25) + 
  geom_segment(aes(x = x - bar_width * 0.1, xend = x + bar_width * 0.1, y = mean + sd, yend = mean + sd), 
               data = error_df , inherit.aes = FALSE, linewidth = 0.25) 

p2 <- p2 + scale_y_continuous(
  limits = c(0, NA),
  expand = expansion(mult = c(0, 0.1))
)

# Combine Figure 2i and Figure 2j panels.
combined_p <- wrap_plots(
  list(
    "inter_group" = p2 + 
      labs(tag = "a") +
      theme(
        legend.key.width = unit(15, "pt"),
        legend.key.height = unit(5, "pt"),
        legend.key.spacing.y = unit(2,"pt"),
        legend.title.position = "left",
        legend.position = "top",
        legend.box.spacing = unit(3, "pt"), 
        plot.tag = element_text(size = 10, margin = ggplot2::margin(0,0,0,0), face = "bold"),  
        plot.margin = ggplot2::margin(0,0,0,0)
      ),
    "intra_group" = p1 + 
      labs(tag = "b") +
      theme(
        legend.key.width = unit(12, "pt"),
        legend.key.height = unit(5, "pt"),
        legend.key.spacing.y = unit(2,"pt"),
        legend.title.position = "left",
        legend.position = "top",
        legend.box.spacing = unit(3, "pt"), 
        plot.tag = element_text(size = 10, margin = ggplot2::margin(0,0,0,0), face = "bold"), 
        plot.margin = ggplot2::margin(0,0,0,0)
      )
  ), 
  ncol = 1,
) + 
  plot_layout(axes = "collect")
ggsave(file.path(COMPARISON_FIGURE_DIR, "angle_diff_test_by_ivdtype.png"), combined_p, width = 180, height = 180*9/16, units = "mm", dpi = 300)




