library(ggplot2)
library(patchwork)
library(scales)
library(latex2exp)
library(reshape2)


options(ggplot2.discrete.fill = NULL)
options(ggplot2.discrete.colour = NULL)


fdr_emp <- read.csv("FDR_results_emp.csv", row.names = 1)
fdr_drm <- read.csv("FDR_results_drm.csv", row.names = 1)
error_emp_normal <- read.csv("error_emp_normal.csv", row.names = 1)
error_drm_normal <- read.csv("error_drm_normal.csv", row.names = 1)
error_emp_gamma <- read.csv("error_emp_gamma.csv", row.names = 1)
error_drm_gamma <- read.csv("error_drm_gamma.csv", row.names = 1)
error_emp_twocompnormal <- read.csv("error_emp_twocompnormal.csv", row.names = 1)
error_drm_twocompnormal <- read.csv("error_drm_twocompnormal.csv", row.names = 1)


calculate_limits <- function(values) {
  max_val <- max(values, na.rm = TRUE)
  if (max_val > 0.4) {
    upper_limit <- ceiling(max_val * 10) / 10
    breaks <- seq(0, upper_limit, by = 0.1)
  } else {
    upper_limit <- ceiling(max_val * 20) / 20
    breaks <- seq(0, upper_limit, by = 0.05)
  }
  return(list(limits = c(0, upper_limit), breaks = breaks))
}

plot_error_and_fdr <- function(error_emp, error_drm, fdr_emp, fdr_drm, dist_name, fdr_colname) {
  df <- data.frame(
    n = factor(as.character(rownames(error_emp)), levels = as.character(rownames(error_emp))),
    Type_I_emp = error_emp[,1],
    Type_II_emp = error_emp[,2],
    Type_I_drm = error_drm[,1],
    Type_II_drm = error_drm[,2],
    FDR_emp = fdr_emp[[fdr_colname]],
    FDR_drm = fdr_drm[[fdr_colname]]
  )
  
  y1 <- calculate_limits(c(df$Type_I_emp, df$Type_I_drm, 0.05))
  y2 <- calculate_limits(c(df$Type_II_emp, df$Type_II_drm))
  y3 <- calculate_limits(c(df$FDR_emp, df$FDR_drm, 0.05))
  
  box_theme <- theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(0.15, "cm"),
      axis.text.x = element_text(margin = margin(t = 3)),
      axis.text.y = element_text(margin = margin(r = 3)),
      plot.title = element_text(hjust = 0.5, size = 12),
      legend.position = "bottom"
    )
  
  method_colors <- c("EMP" = "#1b365d", "DRM(ours)" = "#990000")
  method_linetypes <- c("EMP" = "dashed", "DRM(ours)" = "solid")
  method_shapes <- c("EMP" = 21, "DRM(ours)" = 16)
  
  # Type I error
  df_type1 <- melt(df[, c("n", "Type_I_emp", "Type_I_drm")], id.vars = "n",
                   variable.name = "Method", value.name = "Error")
  df_type1$Method <- factor(df_type1$Method, levels = c("Type_I_emp", "Type_I_drm"),
                            labels = c("EMP", "DRM(ours)"))
  
  p1 <- ggplot(df_type1, aes(x = n, y = Error, 
                             color = Method, linetype = Method, shape = Method, group = Method)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2, fill = "white") +
    geom_hline(yintercept = 0.05, linetype = "dashed") +
    scale_y_continuous(limits = y1$limits, breaks = y1$breaks) +
    scale_color_manual(values = method_colors) +
    scale_linetype_manual(values = method_linetypes) +
    scale_shape_manual(values = method_shapes) +
    box_theme +
    labs(title = TeX(paste0(dist_name, " - Rejection Rate under True $H_0$")),
         x = "Sample size n", y = "Rejection Rate", 
         color = "Method", linetype = "Method", shape = "Method")
  
  # Type II error
  df_type2 <- melt(df[, c("n", "Type_II_emp", "Type_II_drm")], id.vars = "n",
                   variable.name = "Method", value.name = "Error")
  df_type2$Method <- factor(df_type2$Method, levels = c("Type_II_emp", "Type_II_drm"),
                            labels = c("EMP", "DRM(ours)"))
  
  p2 <- ggplot(df_type2, aes(x = n, y = Error, 
                             color = Method, linetype = Method, shape = Method, group = Method)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2, fill = "white") +
    scale_y_continuous(limits = y2$limits, breaks = y2$breaks) +
    scale_color_manual(values = method_colors) +
    scale_linetype_manual(values = method_linetypes) +
    scale_shape_manual(values = method_shapes) +
    box_theme +
    labs(title = TeX(paste0(dist_name, " - Rejection Rate under False $H_0$")),
         x = "Sample size n", y = "Rejection Rate", 
         color = "Method", linetype = "Method", shape = "Method")
  
  # FDR
  df_fdr <- melt(df[, c("n", "FDR_emp", "FDR_drm")], id.vars = "n",
                 variable.name = "Method", value.name = "FDR")
  df_fdr$Method <- factor(df_fdr$Method, levels = c("FDR_emp", "FDR_drm"),
                          labels = c("EMP", "DRM(ours)"))
  
  p3 <- ggplot(df_fdr, aes(x = n, y = FDR, 
                           color = Method, linetype = Method, shape = Method, group = Method)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2, fill = "white") +
    geom_hline(yintercept = 0.05, linetype = "dashed") +
    scale_y_continuous(limits = y3$limits, breaks = y3$breaks) +
    scale_color_manual(values = method_colors) +
    scale_linetype_manual(values = method_linetypes) +
    scale_shape_manual(values = method_shapes) +
    box_theme +
    labs(title = paste(dist_name, "- FDR"),
         x = "Sample size n", y = "False Discovery Rate", 
         color = "Method", linetype = "Method", shape = "Method")
  
  combined_plot <- (p1 + p2 + p3 + plot_layout(ncol = 3)) &
    theme(legend.position = "bottom", legend.key.width = unit(2, "cm"))
  
  return(combined_plot)
}


plot_normal <- plot_error_and_fdr(error_emp_normal, error_drm_normal, fdr_emp, fdr_drm, "Normal", "FDR_normal")
plot_gamma <- plot_error_and_fdr(error_emp_gamma, error_drm_gamma, fdr_emp, fdr_drm, "Gamma", "FDR_gamma")
plot_twocompnormal <- plot_error_and_fdr(error_emp_twocompnormal, error_drm_twocompnormal, fdr_emp, fdr_drm, "Two-Component Normal", "FDR_twocompnormal")

ggsave("errors_Normal.png", plot = plot_normal, width = 15, height = 5, dpi = 300)
ggsave("errors_Gamma.png", plot = plot_gamma, width = 15, height = 5, dpi = 300)
ggsave("errors_Twocompnormal.png", plot = plot_twocompnormal, width = 15, height = 5, dpi = 300)

cat("All figures have been saved！\n")

