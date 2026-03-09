n_list <- c(50, 100, 200, 500, 1000)

# density ratio model
error_twocompnormal <- matrix(0, nrow = 5, ncol = 2)
error_gamma <- matrix(0, nrow = 5, ncol = 2)
error_normal <- matrix(0, nrow = 5, ncol = 2)

for (i in 1:5) {
  n <- n_list[i]
  
  p_twocompnormal <- read.csv(paste0("p_value_drm_twocompnormal_", n, ".csv"))[[1]]
  type1_error_two <- mean(p_twocompnormal[1:1000] <= 0.05)
  type2_error_two <- mean(p_twocompnormal[1001:2000] <= 0.05)
  error_twocompnormal[i, 1] <- type1_error_two
  error_twocompnormal[i, 2] <- type2_error_two
  
  p_gamma <- read.csv(paste0("p_value_drm_gamma_", n, ".csv"))[[1]]
  type1_error_gamma <- mean(p_gamma[1:1000] <= 0.05)
  type2_error_gamma <- mean(p_gamma[1001:2000] <= 0.05)
  error_gamma[i, 1] <- type1_error_gamma
  error_gamma[i, 2] <- type2_error_gamma
  
  p_normal <- read.csv(paste0("p_value_drm_normal_", n, ".csv"))[[1]]
  type1_error_normal <- mean(p_normal[1:1000] <= 0.05)
  type2_error_normal <- mean(p_normal[1001:2000] <= 0.05)
  error_normal[i, 1] <- type1_error_normal
  error_normal[i, 2] <- type2_error_normal
}

colnames(error_twocompnormal) <- c("rejection_rate_TN", "rejection_rate_FN")
rownames(error_twocompnormal) <- n_list

colnames(error_gamma) <- c("rejection_rate_TN", "rejection_rate_FN")
rownames(error_gamma) <- n_list

colnames(error_normal) <- c("rejection_rate_TN", "rejection_rate_FN")
rownames(error_normal) <- n_list

# save files
write.csv(error_twocompnormal, "error_drm_twocompnormal.csv", row.names = TRUE)
write.csv(error_gamma, "error_drm_gamma.csv", row.names = TRUE)
write.csv(error_normal, "error_drm_normal.csv", row.names = TRUE)

cat("All results saved：error_drm_twocompnormal.csv、error_drm_gamma.csv、error_drm_normal.csv\n")


# empirical method
error_twocompnormal <- matrix(0, nrow = 5, ncol = 2)
error_gamma <- matrix(0, nrow = 5, ncol = 2)
error_normal <- matrix(0, nrow = 5, ncol = 2)

for (i in 1:5) {
  n <- n_list[i]
  
  p_twocompnormal <- read.csv(paste0("p_value_emp_twocompnormal_", n, ".csv"))[[1]]
  type1_error_two <- mean(p_twocompnormal[1:1000] <= 0.05)
  type2_error_two <- mean(p_twocompnormal[1001:2000] <= 0.05)
  error_twocompnormal[i, 1] <- type1_error_two
  error_twocompnormal[i, 2] <- type2_error_two
  
  p_gamma <- read.csv(paste0("p_value_emp_gamma_", n, ".csv"))[[1]]
  type1_error_gamma <- mean(p_gamma[1:1000] <= 0.05)
  type2_error_gamma <- mean(p_gamma[1001:2000] <= 0.05)
  error_gamma[i, 1] <- type1_error_gamma
  error_gamma[i, 2] <- type2_error_gamma
  
  p_normal <- read.csv(paste0("p_value_emp_normal_", n, ".csv"))[[1]]
  type1_error_normal <- mean(p_normal[1:1000] <= 0.05)
  type2_error_normal <- mean(p_normal[1001:2000] <= 0.05)
  error_normal[i, 1] <- type1_error_normal
  error_normal[i, 2] <- type2_error_normal
}

colnames(error_twocompnormal) <- c("rejection_rate_TN", "rejection_rate_FN")
rownames(error_twocompnormal) <- n_list

colnames(error_gamma) <- c("rejection_rate_TN", "rejection_rate_FN")
rownames(error_gamma) <- n_list

colnames(error_normal) <- c("rejection_rate_TN", "rejection_rate_FN")
rownames(error_normal) <- n_list

# save files
write.csv(error_twocompnormal, "error_emp_twocompnormal.csv", row.names = TRUE)
write.csv(error_gamma, "error_emp_gamma.csv", row.names = TRUE)
write.csv(error_normal, "error_emp_normal.csv", row.names = TRUE)

cat("All results saved：error_emp_twocompnormal.csv、error_emp_gamma.csv、error_emp_normal.csv\n")
