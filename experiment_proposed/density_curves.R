library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)
library(patchwork)

# set file path
file_path <- "results_all.xlsx"

# define column names
cols_mse <- c("SegRNN_96_96_mse", "ModernTCN_96_96_mse",
              "iTransformer_96_96_mse", "MTST_96_96_mse", 
              "DLinear_96_96_mse", "NLinear_96_96_mse", 
              "MTSMixer_96_96_mse")

cols_atpe <- c("SegRNN_96_96_ATPE", "ModernTCN_96_96_ATPE",
               "iTransformer_96_96_ATPE", "MTST_96_96_ATPE", 
               "DLinear_96_96_ATPE", "NLinear_96_96_ATPE", 
               "MTSMixer_96_96_ATPE")

cols_ms <- c("SegRNN_96_96_ms/sample", "ModernTCN_96_96_ms/sample",
             "iTransformer_96_96_ms/sample", "MTST_96_96_ms/sample", 
             "DLinear_96_96_ms/sample", "NLinear_96_96_ms/sample", 
             "MTSMixer_96_96_ms/sample")

model_names <- c("SegRNN", "ModernTCN", "iTransformer", 
                 "MTST", "DLinear", "NLinear", "MTSMixer")

# choose colors
model_colors <- c("#1B9E77",
                  "#D95F02",
                  "#7570B3",
                  "#66A61E",
                  "#E7298A",
                  "#E6AB02",
                  "grey")

# read data
data_all <- read_excel(file_path)

title_settings <- theme(
  plot.title = element_text(size = 16)
)

# 1. MSE density plot
df_mse <- data_all %>%
  select(all_of(cols_mse)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "model_raw",
    values_to = "value"
  ) %>%
  mutate(model = factor(model_raw, levels = cols_mse, labels = model_names),
         metric = "MSE")

p_mse <- ggplot(df_mse, aes(x = value, fill = model)) +
  geom_density(alpha = 0.5, adjust = 1.5) +
  scale_fill_manual(values = model_colors) +
  labs(x = "MSE", y = "Density", title = "Prediction Error (MSE)") +
  theme_minimal() +
  title_settings +
  theme(
    legend.position = "none",  
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title.y = element_text(size = 14) 
  ) +
  scale_x_continuous(trans = "identity")

# 2. ATPE density plot
df_atpe <- data_all %>%
  select(all_of(cols_atpe)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "model_raw",
    values_to = "value"
  ) %>%
  mutate(model = factor(model_raw, levels = cols_atpe, labels = model_names),
         metric = "ATPE")

p_atpe <- ggplot(df_atpe, aes(x = value + 1e-6, fill = model)) +  
  geom_density(alpha = 0.5, adjust = 1.5) +
  scale_fill_manual(values = model_colors) +
  labs(x = "ATPE", y = "", title = "Training Efficiency (ATPE)") +
  theme_minimal() +
  title_settings +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  scale_x_continuous(trans = "log10")  

# 3. AIT density plot
df_ms <- data_all %>%
  select(all_of(cols_ms)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "model_raw",
    values_to = "value"
  ) %>%
  mutate(model = factor(model_raw, levels = cols_ms, labels = model_names),
         metric = "Inference Time")

p_ms <- ggplot(df_ms, aes(x = value + 1e-6, fill = model)) +
  geom_density(alpha = 0.5, adjust = 1.5) +
  scale_fill_manual(values = model_colors) +
  labs(x = "AIT", y = "", 
       title = "Inference Speed (AIT)") +
  theme_minimal() +
  title_settings +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  scale_x_continuous(trans = "log10")  

# extract the shared legend
legend_plot <- ggplot(df_mse, aes(x = value, fill = model)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = model_colors) +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 14))


model_legend <- cowplot::get_legend(legend_plot)

# combine the three figure
combined_plot <- (p_mse | p_atpe | p_ms) / 
  model_legend +
  plot_layout(heights = c(10, 1))

# save the combined figure
ggsave("density_curves.png", plot = combined_plot, 
       width = 16, height = 6, dpi = 300, bg = "white")

