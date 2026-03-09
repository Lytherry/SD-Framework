# Normal distribution
case1_params <- list(
  F1 = structure(c(0, 0.5**2), type = "normal"),
  G1 = structure(c(0, 2**2), type = "normal"),
  F2 = structure(c(0, 0.5**2), type = "normal"),
  G2 = structure(c(0.75, 2**2), type = "normal")
)

# Gamma distribution
case2_params <- list(
  F1 = structure(c(12, 3), type = "gamma"),
  G1 = structure(c(2, 0.432), type = "gamma"),
  F2 = structure(c(12, 3), type = "gamma"),
  G2 = structure(c(2, 0.333), type = "gamma")
)

# Two-Component Normal Mixture distribution
mixnorm_F <- structure(
  list(
    list(mean = 19, sd = 0.5),
    list(mean = 21, sd = 0.5)
  ), 
  type = "mixture_normal"
)

mixnorm_G1 <- structure(
  list(
    list(mean = 17, sd = 2),
    list(mean = 23, sd = 2)
  ), 
  type = "mixture_normal",
  weights = c(0.5, 0.5)
)

mixnorm_G2 <- structure(
  list(
    list(mean = 17, sd = 2),
    list(mean = 23, sd = 2)
  ), 
  type = "mixture_normal",
  weights = c(0.34, 0.66)
)

case3_params <- list(
  F1 = mixnorm_F,
  G1 = mixnorm_G1,
  F2 = mixnorm_F,
  G2 = mixnorm_G2
)


generate_combined_quantile <- function(F_params, G_params, p_seq = seq(0.000000001, 0.999999999, length.out = 1000)) {
  F_type <- attr(F_params, "type")
  G_type <- attr(G_params, "type")
  
  F_quant <- switch(F_type,
                    "normal" = qnorm(p_seq, F_params[1], sqrt(F_params[2])),
                    "gamma" = qgamma(p_seq, shape = F_params[1], rate = F_params[2]),
                    "mixture_normal" = {
                      weights <- attr(F_params, "weights", exact = TRUE)
                      if (is.null(weights)) weights <- c(0.5, 0.5)
                      sapply(p_seq, function(p) {
                        uniroot(function(x) {
                          weights[1] * pnorm(x, F_params[[1]]$mean, F_params[[1]]$sd) +
                            weights[2] * pnorm(x, F_params[[2]]$mean, F_params[[2]]$sd) - p
                        }, interval = c(min(F_params[[1]]$mean, F_params[[2]]$mean) - 6*max(F_params[[1]]$sd,F_params[[2]]$sd),
                                        max(F_params[[1]]$mean, F_params[[2]]$mean) + 6*max(F_params[[1]]$sd,F_params[[2]]$sd)))$root
                      })
                    }
  )
  
  G_quant <- switch(G_type,
                    "normal" = qnorm(p_seq, G_params[1], sqrt(G_params[2])),
                    "gamma" = qgamma(p_seq, shape = G_params[1], rate = G_params[2]),
                    "mixture_normal" = {
                      weights <- attr(G_params, "weights", exact = TRUE)
                      if (is.null(weights)) weights <- c(0.5, 0.5)
                      sapply(p_seq, function(p) {
                        uniroot(function(x) {
                          weights[1] * pnorm(x, G_params[[1]]$mean, G_params[[1]]$sd) +
                            weights[2] * pnorm(x, G_params[[2]]$mean, G_params[[2]]$sd) - p
                        }, interval = c(min(G_params[[1]]$mean, G_params[[2]]$mean) - 6*max(G_params[[1]]$sd,G_params[[2]]$sd),
                                        max(G_params[[1]]$mean, G_params[[2]]$mean) + 6*max(G_params[[1]]$sd,G_params[[2]]$sd)))$root
                      })
                    }
  )
  
  list(
    F = data.frame(p = p_seq, q = F_quant),
    G = data.frame(p = p_seq, q = G_quant)
  )
}


plot_comparison_quantile_base <- function(F1_params, G1_params, F2_params, G2_params, 
                                          case_name, filename, y_range = NULL,
                                          vlines_dim1 = NULL, vlines_dim2 = NULL) {
  dim1_data <- generate_combined_quantile(F1_params, G1_params)
  dim2_data <- generate_combined_quantile(F2_params, G2_params)
  
  if (is.null(y_range)) {
    all_q <- c(dim1_data$F$q, dim1_data$G$q, dim2_data$F$q, dim2_data$G$q)
    y_range <- range(all_q, na.rm = TRUE)
  }
  
  par(mfrow = c(1, 2), 
      mar = c(4.5, 3.5, 3.5, 1), 
      las = 1, 
      tcl = -0.3, 
      mgp = c(2, 0.5, 0), 
      oma = c(0, 0, 0, 0)
  )
  
  # --- Dimension 1 ---
  plot(dim1_data$F$p, dim1_data$F$q, type = "l", col = "#1b365d", lwd = 2, lty=3, 
       xlab = "t", ylab = "Quantile", main = paste(case_name, "- Dimension 1"),
       ylim = y_range)
  lines(dim1_data$G$p, dim1_data$G$q, col = "#990000", lwd = 2)
  legend("topleft", col = c("#1b365d", "#990000"), lty = c(3,1), lwd = 2,
         legend = c(expression(F[1]), expression(G[1])))
  
  if (!is.null(vlines_dim1)) {
    for (p_target in vlines_dim1) {
      q_val <- approx(dim1_data$F$p, dim1_data$F$q, xout = p_target)$y
      y_min <- par("usr")[3]
      segments(p_target, q_val, p_target, y_min, lty = 2, lwd = 2, col = "gray40")
      
      text(x = p_target + 0.02,
           y = (q_val + y_min)/2,
           labels = paste0("t=", round(p_target, 4)),
           pos = 4,
           cex = 1.2, col = "gray20")
    }
  }
  
  # --- Dimension 2 ---
  plot(dim2_data$F$p, dim2_data$F$q, type = "l", col = "#1b365d", lwd = 2, lty=3, 
       xlab = "t", ylab = "Quantile", main = paste(case_name, "- Dimension 2"),
       ylim = y_range)
  lines(dim2_data$G$p, dim2_data$G$q, col = "#990000", lwd = 2)
  legend("topleft", col = c("#1b365d", "#990000"), lty = c(3,1), lwd = 2,
         legend = c(expression(F[2]), expression(G[2])))
  
  if (!is.null(vlines_dim2)) {
    for (p_target in vlines_dim2) {
      q_val <- approx(dim2_data$F$p, dim2_data$F$q, xout = p_target)$y
      y_min <- par("usr")[3]
      segments(p_target, q_val, p_target, y_min, lty = 2, lwd = 2, col = "gray40")
      
      text(x = p_target + 0.02, 
           y = (q_val + y_min)/2, 
           labels = paste0("t=", round(p_target, 4)),
           pos = 4,
           cex = 1.2, col = "gray20")
    }
  }
  
  dev.copy(png, filename = filename, width = 8.5, height = 5, units = "in", res = 300)
  dev.off()
  par(mfrow = c(1, 1))
}


plot_comparison_quantile_base(case1_params$F1, case1_params$G1, 
                              case1_params$F2, case1_params$G2, 
                              "Normal", "quantile_normal.png",
                              y_range = c(-2.5,2.5),
                              vlines_dim1 = c(0.5),
                              vlines_dim2 = c(0.3085))

plot_comparison_quantile_base(case2_params$F1, case2_params$G1, 
                              case2_params$F2, case2_params$G2, 
                              "Gamma", "quantile_gamma.png",
                              y_range = c(1,13),
                              vlines_dim1 = c(0.5),
                              vlines_dim2 = c(0.3054))

plot_comparison_quantile_base(case3_params$F1, case3_params$G1, 
                              case3_params$F2, case3_params$G2, 
                              "Two-Component Normal", "quantile_twocompnormal.png",
                              y_range = c(17,23),
                              vlines_dim1 = c(0.5),
                              vlines_dim2 = c(0.3097))

cat("All figures have been saved！\n")

