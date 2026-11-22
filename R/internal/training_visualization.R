# DiscSeek | Parallel V-Net training and validation curves
# Outputs: Supplementary Figure S4.

source("R/00_config.R")
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
# 1. Read V-Net training logs ----
metric_df <- read.csv(file.path(INPUT_DIR, "training", "loss_metric_curves.csv"))
plot_list <- list()

# Figure S4a | Segmentation V-Net training loss ----
plot_list[['train_seg_loss']] <- metric_df %>% 
  filter(Network == 'segmentation', Dataset == 'training') %>%
  group_by(Epoch) %>%
  summarise(Loss = mean(Loss)) %>% 
  filter(Epoch < 137) %>%
  ggplot(aes(x=Epoch, y=Loss)) +
  geom_line(linewidth=0.4, colour="#EC5D3B") +
  scale_x_continuous(breaks = seq(0, 136, by=40)) + 
  labs(x='Epoch', y='Training Focal Loss') + 
  coord_cartesian(expand = FALSE) +
  labs(tag = "a") +
  theme(axis.title = element_text(size = 10), axis.text = element_text(size = 8),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.25),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.25),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 8, margin = ggplot2::margin(t=0,b=1)),
        plot.tag = element_text(size = 10, margin = ggplot2::margin(0,0,0,0), face = "bold"), 
        plot.margin = ggplot2::margin(0,0,0,0)
  )
# Figure S4d | Identification V-Net training loss ----
plot_list[['train_ide_loss']] <- metric_df %>% filter(Network == 'identification', Dataset == 'training') %>%
  group_by(Epoch) %>%
  summarise(Loss = mean(Loss)) %>% 
  filter(Epoch < 104) %>%
  ggplot(aes(x=Epoch, y=Loss)) +
  geom_line(linewidth=0.4, colour="#EC5D3B") +
  scale_x_continuous(breaks = seq(0, 103, by=20)) + 
  labs(x='Epoch', y='Training MSE Loss') + 
  coord_cartesian(expand = FALSE) +
  labs(tag = "d") +
  theme(axis.title = element_text(size = 10), axis.text = element_text(size = 8),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.25),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.25),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 8, margin = ggplot2::margin(t=0,b=1)),
        plot.tag = element_text(size = 10, margin = ggplot2::margin(0,0,0,0), face = "bold"), 
        plot.margin = ggplot2::margin(0,0,0,0)
  )
# Figure S4b | Segmentation V-Net validation loss ----
plot_list[['val_seg_loss']] <- metric_df %>% filter(Network == 'segmentation', Dataset == 'validation') %>%
  filter(Epoch < 137) %>%
  ggplot(aes(x=Epoch, y=Loss)) +
  geom_line(linewidth=0.4, colour="#EC5D3B") +
  scale_x_continuous(breaks = seq(0, 136, by=40)) + 
  labs(x='Epoch', y='Validation Focal Loss') + 
  coord_cartesian(expand = FALSE) +
  labs(tag = "b") +
  theme(axis.title = element_text(size = 10), axis.text = element_text(size = 8),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.25),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.25),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 8, margin = ggplot2::margin(t=0,b=1)),
        plot.tag = element_text(size = 10, margin = ggplot2::margin(0,0,0,0), face = "bold"), 
        plot.margin = ggplot2::margin(0,0,0,0)
  )
# Figure S4e | Identification V-Net validation loss ----
plot_list[['val_ide_loss']] <- metric_df %>% filter(Network == 'identification', Dataset == 'validation') %>%
  filter(Epoch < 104) %>%
  ggplot(aes(x=Epoch, y=Loss)) +
  geom_line(linewidth=0.4, colour="#EC5D3B") +
  scale_x_continuous(breaks = seq(0, 103, by=20)) + 
  labs(x='Epoch', y='Validation MSE Loss') + 
  coord_cartesian(expand = FALSE) +
  labs(tag = "e") +
  theme(axis.title = element_text(size = 10), axis.text = element_text(size = 8),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.25),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.25),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 8, margin = ggplot2::margin(t=0,b=1)),
        plot.tag = element_text(size = 10, margin = ggplot2::margin(0,0,0,0), face = "bold"), 
        plot.margin = ggplot2::margin(0,0,0,0)
  )
# Figure S4c | Segmentation V-Net validation Dice ----
plot_list[['val_seg_dice']] <- metric_df %>% filter(Network == 'segmentation', Dataset == 'validation') %>%
  filter(Epoch < 137) %>%
  ggplot(aes(x=Epoch, y=Dice_ID_Rate)) +
  geom_line(linewidth=0.4, colour="#EC5D3B") +
  scale_x_continuous(breaks = seq(0, 136, by=40)) + 
  labs(x='Epoch', y='Validation Dice') + 
  coord_cartesian(expand = FALSE) +
  labs(tag = "c") +
  theme(axis.title = element_text(size = 10), axis.text = element_text(size = 8),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.25),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.25),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 8, margin = ggplot2::margin(t=0,b=1)),
        plot.tag = element_text(size = 10, margin = ggplot2::margin(0,0,0,0), face = "bold"), 
        plot.margin = ggplot2::margin(0,0,0,0)
  )
# Figure S4f | Identification V-Net validation identified rate ----
plot_list[['val_ide_rate']] <- metric_df %>% filter(Network == 'identification', Dataset == 'validation') %>%
  filter(Epoch < 104) %>%
  ggplot(aes(x=Epoch, y=Dice_ID_Rate)) +
  geom_line(linewidth=0.4, colour="#EC5D3B") +
  scale_x_continuous(breaks = seq(0, 103, by=20)) + 
  labs(x='Epoch', y='Validation Identified Rate') + 
  coord_cartesian(expand = FALSE) +
  labs(tag = "f") +
  theme(axis.title = element_text(size = 10), axis.text = element_text(size = 8),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.25),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.25),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 8, margin = ggplot2::margin(t=0,b=1)),
        plot.tag = element_text(size = 10, margin = ggplot2::margin(0,0,0,0), face = "bold"), 
        plot.margin = ggplot2::margin(0,0,0,0)
  )
# Combine Supplementary Figure S4 panels ----
combined_p <- wrap_plots(
  plot_list, 
  ncol = 2,
) + 
  plot_layout(axes = "collect_x")

ggsave(file.path(COMPARISON_FIGURE_DIR, "network_curves.eps"), combined_p, width = 180, height = 180*12/16, units = "mm", dpi = 300)
