n_list <- c(50, 100, 200, 500, 1000)

fdr_twocompnormal <- numeric(5)
fdr_gamma <- numeric(5)
fdr_normal <- numeric(5)

compute_fdr <- function(p_values, threshold = 0.05) {
  reject <- which(p_values <= threshold)
  if (length(reject) == 0) {
    return(0)
  }
  false_positives <- sum(reject <= 1000)
  fdr <- false_positives / length(reject)
  return(fdr)
}

for (i in 1:5) {
  n <- n_list[i]
  
  adj_p_twocompnormal <- read.csv(paste0("adjusted_p_value_drm_twocompnormal_", n, ".csv"))[[1]]
  fdr_twocompnormal[i] <- compute_fdr(adj_p_twocompnormal)
  
  adj_p_gamma <- read.csv(paste0("adjusted_p_value_drm_gamma_", n, ".csv"))[[1]]
  fdr_gamma[i] <- compute_fdr(adj_p_gamma)
  
  adj_p_normal <- read.csv(paste0("adjusted_p_value_drm_normal_", n, ".csv"))[[1]]
  fdr_normal[i] <- compute_fdr(adj_p_normal)
}

fdr_results <- data.frame(
  n = n_list,
  FDR_twocompnormal = fdr_twocompnormal,
  FDR_gamma = fdr_gamma,
  FDR_normal = fdr_normal
)

write.csv(fdr_results, "FDR_results_drm.csv", row.names = FALSE)

cat("FDR calculation completed, and the results have been saved to FDR_results_drm.csv\n")


for (i in 1:5) {
  n <- n_list[i]
  
  adj_p_twocompnormal <- read.csv(paste0("adjusted_p_value_emp_twocompnormal_", n, ".csv"))[[1]]
  fdr_twocompnormal[i] <- compute_fdr(adj_p_twocompnormal)
  
  adj_p_gamma <- read.csv(paste0("adjusted_p_value_emp_gamma_", n, ".csv"))[[1]]
  fdr_gamma[i] <- compute_fdr(adj_p_gamma)
  
  adj_p_normal <- read.csv(paste0("adjusted_p_value_emp_normal_", n, ".csv"))[[1]]
  fdr_normal[i] <- compute_fdr(adj_p_normal)
}

fdr_results <- data.frame(
  n = n_list,
  FDR_twocompnormal = fdr_twocompnormal,
  FDR_gamma = fdr_gamma,
  FDR_normal = fdr_normal
)

write.csv(fdr_results, "FDR_results_emp.csv", row.names = FALSE)

cat("FDR calculation completed, and the results have been saved to FDR_results_emp.csv\n")
