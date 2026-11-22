# DiscSeek | Internal retrospective technical validation
# Outputs: Table 1, Supplementary Figures S2–S3.

# 1. Load internal validation data ----
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

origin_internal_data <- read_csv(file.path(INPUT_DIR, "internal", "internal_planning_error.csv"))
# Input variable definitions:
# ID: unique internal test examination identifier.
# SAG_dirET, SAG_cenEA, SAG_dirEA: sagittal inter-technologist angle difference, AI-to-expert FoV center offset, and AI-to-expert angle offset, respectively.
# COR_dirET, COR_cenEA, COR_dirEA: corresponding coronal measurements.
# AXI_dirET, AXI_cenEA, AXI_dirEA: corresponding axial measurements.
# SAG_DISC: visible disc count in the best-viewed sagittal image.
# COR_DISC: visible disc count in the best-viewed coronal image.
# AXI_DISC: visible disc count in the scout volume.
# Qualified.Mask: technologist-assessed segmentation mask qualification.
# Quality.Factor: score calculated using the predefined mask qualification criteria.

# 2. Supplementary Figure S2 | Visible disc coverage ----
# Figure S2a | Compare sagittal best-viewed images and scouts.
dis_table <- origin_internal_data %>% 
  rowwise() %>% 
  mutate(SAG = min(SAG_DISC, AXI_DISC), 
         COR_DISC = min(COR_DISC, AXI_DISC)) %>%  
  mutate(SAG_DISC = min(SAG_DISC, 6), 
         COR_DISC = min(COR_DISC, 6),
         AXI_DISC = min(AXI_DISC, 6)) %>% 
  group_by(SAG_DISC, AXI_DISC) %>% 
  summarise(count = n(), .groups = "drop")
p1 <- dis_table %>%
  complete(SAG_DISC = unique(SAG_DISC), AXI_DISC = unique(AXI_DISC), fill = list(count = 0)) %>%
  ggplot(aes(x = SAG_DISC, y = AXI_DISC)) + 
  geom_tile(aes(fill = count), color = "black", show.legend = TRUE) + 
  geom_text(aes(label = count), color = "black", size.unit = "pt", size = 10) + 
  scale_x_continuous(breaks = seq(4, 6, by=1)) + 
  scale_y_continuous(breaks = seq(4, 6, by=1)) + 
  scale_fill_gradientn(
    name = "Image\nCount",
    colors = colors[3:7],
    na.value = "grey99"
  ) + 
  labs(y = "Discs (MR Scout)", x = "Discs (Best Sagittal)") + 
  coord_cartesian(expand = FALSE) + 
  # coord_fixed() + 
  theme(
    plot.title = element_text(size = 8, hjust = 0.5, margin = ggplot2::margin(t=5, b=3, unit = "pt")),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 8, margin = ggplot2::margin(t=0,b=0)),
    legend.title = element_text(size = 10, margin = ggplot2::margin(r=10, unit = "pt")),  # , angle = 90, hjust = 0.5),
    legend.position = "bottom",
    legend.box.spacing = unit(3, "pt"),
    legend.margin = ggplot2::margin(0,0,0,0),
    legend.justification = "center",
    legend.box.just = "center",
    axis.text = element_text(color = "black", size = 8),
    axis.title = element_text(size = 10),
    panel.grid.minor = element_line(colour = "grey", linewidth=0.25),
    panel.grid.major = element_line(colour = "grey", linewidth = 0.25),
    panel.background = element_rect(fill = "white", colour = "grey", linewidth = 0.5),
  )

# Figure S2b | Compare coronal best-viewed images and scouts.
dis_table <- origin_internal_data %>% 
  rowwise() %>% 
  mutate(SAG_DISC = min(SAG_DISC, AXI_DISC), 
         COR_DISC = min(COR_DISC, AXI_DISC)) %>%  
  mutate(SAG_DISC = min(SAG_DISC, 6), 
         COR_DISC = min(COR_DISC, 6),
         AXI_DISC = min(AXI_DISC, 6)) %>% 
  group_by(COR_DISC, AXI_DISC) %>% 
  summarise(count = n(), .groups = "drop")
p2 <- dis_table %>%
  complete(COR_DISC = unique(COR_DISC), AXI_DISC = unique(AXI_DISC), fill = list(count = 0)) %>%
  ggplot(aes(x = COR_DISC, y = AXI_DISC)) + 
  geom_tile(aes(fill = count), color = "black", show.legend = TRUE) + 
  geom_text(aes(label = count), color = "black", size.unit = "pt", size = 10) + 
  scale_x_continuous(breaks = seq(4, 6, by=1)) + 
  scale_y_continuous(breaks = seq(4, 6, by=1)) + 
  scale_fill_gradientn(
    name = "Image\nCount",
    colors = colors[3:7],
    na.value = "grey99"
  ) + 
  labs(y = "Discs (MR Scout)", x = "Discs (Best Coronal)") + 
  coord_cartesian(expand = FALSE) + 
  # coord_fixed() + 
  theme(
    plot.title = element_text(size = 8, hjust = 0.5, margin = ggplot2::margin(t=5, b=3, unit = "pt")),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 8, margin = ggplot2::margin(t=0,b=0)),
    legend.title = element_text(size = 10, margin = ggplot2::margin(r=10, unit = "pt")),  # , angle = 90, hjust = 0.5),
    legend.position = "bottom",
    legend.box.spacing = unit(3, "pt"),
    legend.margin = ggplot2::margin(0,0,0,0),
    legend.justification = "center",
    legend.box.just = "center",
    axis.text = element_text(color = "black", size = 8),
    axis.title = element_text(size = 10),
    panel.grid.minor = element_line(colour = "grey", linewidth=0.25),
    panel.grid.major = element_line(colour = "grey", linewidth = 0.25),
    panel.background = element_rect(fill = "white", colour = "grey", linewidth = 0.5),
  )

# Figure S2c | Compare sagittal and coronal best-viewed images.
dis_table <- origin_internal_data %>% 
  rowwise() %>% 
  mutate(SAG_DISC = min(SAG_DISC, AXI_DISC), 
         COR_DISC = min(COR_DISC, AXI_DISC)) %>%  
  mutate(SAG_DISC = min(SAG_DISC, 6), 
         COR_DISC = min(COR_DISC, 6),
         AXI_DISC = min(AXI_DISC, 6)) %>% 
  group_by(SAG_DISC, COR_DISC) %>% 
  summarise(count = n(), .groups = "drop")
p3 <- dis_table %>%
  complete(SAG_DISC = unique(SAG_DISC), COR_DISC = unique(COR_DISC), fill = list(count = 0)) %>%
  ggplot(aes(x = SAG_DISC, y = COR_DISC)) + 
  geom_tile(aes(fill = count), color = "black", show.legend = TRUE) + 
  geom_text(aes(label = count), color = "black", size.unit = "pt", size = 10) + 
  scale_x_continuous(breaks = seq(4, 6, by=1)) + 
  scale_y_continuous(breaks = seq(4, 6, by=1)) + 
  scale_fill_gradientn(
    name = "Image\nCount",
    colors = colors[3:7],
    na.value = "grey99"
  ) + 
  labs(y = "Discs (Best Coronal)", x = "Discs (Best Sagittal)") + 
  coord_cartesian(expand = FALSE) + 
  # coord_fixed() + 
  theme(
    plot.title = element_text(size = 8, hjust = 0.5, margin = ggplot2::margin(t=5, b=3, unit = "pt")),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 8, margin = ggplot2::margin(t=0,b=0)),
    legend.title = element_text(size = 10, margin = ggplot2::margin(r=10, unit = "pt")),  # , angle = 90, hjust = 0.5),
    legend.position = "bottom",
    legend.box.spacing = unit(3, "pt"),
    legend.margin = ggplot2::margin(0,0,0,0),
    legend.justification = "center",
    legend.box.just = "center",
    axis.text = element_text(color = "black", size = 8),
    axis.title = element_text(size = 10),
    panel.grid.minor = element_line(colour = "grey", linewidth=0.25),
    panel.grid.major = element_line(colour = "grey", linewidth = 0.25),
    panel.background = element_rect(fill = "white", colour = "grey", linewidth = 0.5),
  )

# Combine Supplementary Figure S2 panels.
combined_p <- wrap_plots(
  list(
    "p1" = p1 + 
      labs(tag = "b") +
      theme(
        plot.tag = element_text(size = 10, margin = ggplot2::margin(0,0,0,0), face = "bold"),  
        plot.margin = ggplot2::margin(0,0,0,0)
      ),
    "p2" = p2 + 
      labs(tag = "c") +
      theme( 
        plot.tag = element_text(size = 10, margin = ggplot2::margin(0,0,0,0), face = "bold"), 
        plot.margin = ggplot2::margin(0,0,0,0)
      ),
    "p3" = p3 + 
      labs(tag = "d") +
      theme( 
        plot.tag = element_text(size = 10, margin = ggplot2::margin(0,0,0,0), face = "bold"), 
        plot.margin = ggplot2::margin(0,0,0,0)
      )
  ), 
  nrow = 1
) + 
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave(file.path(COMPARISON_FIGURE_DIR, "internal_disc_coverage_diff_test.eps"), combined_p, width = 180, height = 180*5/16, units = "mm", dpi = 300)


# 3. Table 1 and Supplementary Figure S3 | Planning errors ----
internal_data <- origin_internal_data %>% dplyr::select(-Qualified.Mask, -Quality.Factor) %>%
  pivot_longer(setdiff(names(origin_internal_data), c("ID", "Qualified.Mask", "Quality.Factor")), names_to = "Parameter", values_to = "Value") %>% 
  separate(Parameter, into = c("Orientation", "Parameter"), sep = "_")

# Table 1 | Summarise angle and center offsets. ----
internal_data %>% 
  filter(Parameter != "DISC") %>% 
  group_by(Orientation, Parameter) %>% 
  summarise(mean = mean(Value), sd=sd(Value)) %>% 
  mutate(error = sprintf(
    "%.2f (%.2f)",
    mean,
    sd
  )) %>%
  dplyr::select(-c(mean, sd)) %>%
  pivot_wider(names_from = Orientation, values_from = error)

# Figure S3a | Visualize scan angle offset distributions. ----
angle_p <- internal_data %>% 
  filter(str_detect(Parameter, "dir")) %>% 
  mutate(Parameter = ifelse(Parameter == "dirEA", "DiscSeek", "Inter-technologist"), 
         Orientation = case_when(Orientation == "SAG" ~ "Sagittal", 
                                 Orientation == "COR" ~ "Coronal",
                                 Orientation == "AXI" ~ "Axial")) %>%
  mutate(Parameter = factor(Parameter, levels = c("Inter-technologist", "DiscSeek")),
         Orientation = factor(Orientation, levels = c("Coronal", "Sagittal", "Axial"))) %>%
  ggplot(aes(x = Orientation, y = Value, color = Parameter)) + 
  geom_violin(position = position_dodge(1.0), alpha=1.0, width = 1.0, linewidth = 0.3) + 
  geom_jitter(stroke = 0.3, position = position_jitterdodge(jitter.width = 0.4, dodge.width = 1.0), 
              size = 1, alpha = 0.5, show.legend = FALSE) + 
  geom_boxplot(width=0.15, position = position_dodge(1.0), outliers=FALSE, alpha = 0.7, linewidth = 0.3) + 
  labs(x = "Imaging Orientation", y = "Scan Angle Offset [°]", color = NULL) + 
  scale_color_manual(values = color_palettes[c(1,6)]) + 
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
# Compare DiscSeek offsets with inter-technologist differences by orientation.
max_y <- internal_data %>% 
  filter(str_detect(Parameter, "dir")) %>% 
  mutate(Parameter = ifelse(Parameter == "dirEA", "DiscSeek", "Inter-technologist"), 
         Orientation = case_when(Orientation == "SAG" ~ "Sagittal", 
                                 Orientation == "COR" ~ "Coronal",
                                 Orientation == "AXI" ~ "Axial")) %>%
  mutate(Parameter = factor(Parameter, levels = c("Inter-technologist", "DiscSeek")),
         Orientation = factor(Orientation, levels = c("Coronal", "Sagittal", "Axial"))) %>%
  group_by(Orientation) %>% summarise(y_max = max(Value, na.rm = TRUE))

posthoc_df <- internal_data %>% distinct(Orientation) %>%
  pmap_dfr(~{
    orientation <- ..1
    result <- compare_means(Value~Parameter, internal_data %>% filter(str_detect(Parameter, "dir"), Orientation == orientation), 
                            method = "wilcox.test", paired = TRUE) %>% 
      mutate(Orientation = orientation, group1 = "DiscSeek", group2 = "Inter-technologist")
  }) %>% 
  mutate(Orientation = case_when(Orientation == "SAG" ~ "Sagittal", 
                                 Orientation == "COR" ~ "Coronal",
                                 Orientation == "AXI" ~ "Axial")) %>%
  mutate(Orientation = factor(Orientation, levels = c("Coronal", "Sagittal", "Axial")))

posthoc_df <- posthoc_df %>% left_join(max_y, by = c("Orientation")) %>% 
  mutate(x = match(Orientation, levels(posthoc_df$Orientation)), 
         xmin = x-0.25,
         xmax = x+0.25,
         y = y_max * 1.02,
         x_text = (xmin + xmax)/2)
angle_p <- angle_p +
  geom_segment(aes(x = xmin, xend = xmax, y = y, yend = y), 
               data = posthoc_df %>% filter(p.signif != "ns"), inherit.aes = FALSE, linewidth = 0.25) + 
  geom_segment(aes(x = xmin, xend = xmin, y = 0.99 * y, yend = y), 
               data = posthoc_df %>% filter(p.signif != "ns"), inherit.aes = FALSE, linewidth = 0.25) + 
  geom_segment(aes(x = xmax, xend = xmax, y = 0.99 * y, yend = y), 
               data = posthoc_df %>% filter(p.signif != "ns"), inherit.aes = FALSE, linewidth = 0.25) + 
  geom_text(aes(label = p.signif, x = x_text, y = y),
            data = posthoc_df %>% filter(p.signif != "ns"), inherit.aes = FALSE, size = 3, vjust = 0, nudge_y = -0.05)  + 
  scale_y_continuous(limits = c(0, 6.8), expand = expansion(mult = c(0.01,0)))

# Figure S3b | Visualize FoV center offset distributions.
dist_p <- internal_data %>% 
  filter(str_detect(Parameter, "cen")) %>% 
  mutate(Parameter = "DiscSeek", 
         Orientation = case_when(Orientation == "SAG" ~ "Sagittal", 
                                 Orientation == "COR" ~ "Coronal",
                                 Orientation == "AXI" ~ "Axial")) %>%
  mutate(Orientation = factor(Orientation, levels = c("Coronal", "Sagittal", "Axial"))) %>%
  ggplot(aes(x = Orientation, y = Value, color = Orientation)) + 
  geom_violin(position = position_dodge(1.0), alpha=1.0, width = 1.0, show.legend = FALSE, linewidth = 0.3) + 
  geom_jitter(stroke=0.3, position = position_jitterdodge(jitter.width = 1.0, dodge.width = 1.0), 
              size = 1, alpha = 0.5, show.legend = FALSE) + 
  geom_boxplot(width=0.15, position = position_dodge(1.0), outliers=FALSE, alpha = 0.7, show.legend = FALSE, linewidth = 0.3) + 
  labs(x = "Imaging Orientation", y = "FoV Center Offset [mm]") + 
  scale_color_manual(values = color_palettes[c(1,5,6)]) + 
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
        strip.text = element_text(size = 8, margin = ggplot2::margin(t=0,b=1)))  + 
  scale_y_continuous(limits = c(NA, 8), expand = expansion(mult = c(0.01,0)))

# Combine Supplementary Figure S3 panels.
combined_p <- wrap_plots(
  list(
    "angle" = angle_p + 
      labs(tag = "a") +
      theme(
        legend.position = c(0.65,0.95),
        legend.box.spacing = unit(3, "pt"), 
        legend.direction = "horizontal",
        plot.tag = element_text(size = 10, margin = ggplot2::margin(0,0,0,0), face = "bold"),  
        plot.margin = ggplot2::margin(0,0,0,0)
      ),
    "dist" = dist_p + 
      labs(tag = "b") +
      theme( 
        plot.tag = element_text(size = 10, margin = ggplot2::margin(0,0,0,0), face = "bold"), 
        plot.margin = ggplot2::margin(0,0,0,0)
      )
  ), 
  nrow = 1,
  widths = (c(1.5,1.0))
) + 
  plot_layout(axes = "collect")
ggsave(file.path(COMPARISON_FIGURE_DIR, "internal_angle_diff_test.png"), combined_p, width = 180, height = 180*7.5/16, units = "mm", dpi = 300)


# Optional visualization | Segmentation quality factor by mask qualification ----
# qualified_p <- internal_data %>% 
#   mutate(Qualified.Mask = ifelse(Qualified.Mask == "qualified", "Qualified", "Unqualified")) %>%
#   mutate(Qualified.Mask = factor(Qualified.Mask, levels = c("Unqualified", "Qualified"))) %>%
#   ggplot(aes(x = Qualified.Mask, y = Quality.Factor, color = Qualified.Mask)) + 
#   geom_violin(position = position_dodge(1.0), alpha=1.0, width = 1.0, show.legend = FALSE, linewidth = 0.3) + 
#   geom_jitter(stroke=0.3, position = position_jitterdodge(jitter.width = 0.9, dodge.width = 1.0), 
#               size = 1, alpha = 0.5, show.legend = FALSE) + 
#   geom_boxplot(width=0.15, position = position_dodge(1.0), outliers=FALSE, alpha = 0.7, show.legend = FALSE, linewidth = 0.3) + 
#   geom_hline(yintercept  = 0.6, colour = colors[5], linetype = "dashed") + 
#   labs(x = "DiscSeek Segmentation\nRating", y = "Quality Factor") + 
#   scale_color_manual(values = colors[c(1,6)]) + 
#   theme(legend.position = "top", legend.margin = ggplot2::margin(b=0),
#         legend.title = element_text(size = 10), legend.key.height = unit(5, "pt"),
#         legend.spacing = unit(3,"pt"), legend.text = element_text(size = 8), legend.box.spacing = unit(0,"mm"),
#         legend.key = element_rect(linewidth = 0, color = "white"), axis.text = element_text(color = "black", size = 8),
#         axis.title = element_text(size = 10), axis.text.y = element_text(size = 6),
#         panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
#         panel.grid.major = element_line(colour = "grey", linewidth = 0.25),
#         panel.grid.minor = element_line(colour = "grey", linewidth = 0.25),
#         panel.grid.major.x = element_blank(),
#         strip.background = element_rect(fill = "white"),
#         strip.text = element_text(size = 8, margin = ggplot2::margin(t=0,b=1)))
# 
# ggsave(file.path(COMPARISON_FIGURE_DIR, "internal_qualified_diff_test.eps"), qualified_p, width = 180*0.4, height = 180*5/16, units = "mm", dpi = 300)
# 
# Figure S2 | See the manuscript mapping for this analysis block.
# ggsave(file.path(COMPARISON_FIGURE_DIR, "internal_qualified_diff_test_1.png"),
#        wrap_plots(
#          list(
#            "p1" = qualified_p + 
#              labs(tag = "a") +
#              theme(
#                plot.tag = element_text(size = 10, margin = ggplot2::margin(0,0,0,0), face = "bold"),  
#                plot.margin = ggplot2::margin(0,0,0,0)
#              ),
#            "p2" = p1 + 
#              labs(tag = "b") +
#              theme( 
#                plot.tag = element_text(size = 10, margin = ggplot2::margin(0,0,0,0), face = "bold"), 
#                plot.margin = ggplot2::margin(0,0,0,0)
#              ),
#            "p3" = p2 + 
#              labs(tag = "c") +
#              theme( 
#                plot.tag = element_text(size = 10, margin = ggplot2::margin(0,0,0,0), face = "bold"), 
#                plot.margin = ggplot2::margin(0,0,0,0)
#              ),
#            "p4" = p3 + 
#              labs(tag = "d") +
#              theme( 
#                plot.tag = element_text(size = 10, margin = ggplot2::margin(0,0,0,0), face = "bold"), 
#                plot.margin = ggplot2::margin(0,0,0,0)
#              )
#          ), 
#          nrow = 1
#        ) + 
#          plot_layout(guides = "collect") &
#          theme(legend.position = "bottom", legend.key.height = unit(8, units = "pt")), 
#        width = 180, height = 180*5.8/16, units = "mm", dpi = 300)
