library(parallel)
library(readxl)
source('Functions.R')

quantiles <- seq(from=0.001, to=0.999, length.out=999)
BasisFun <- c('x','x^2','log(x)')

B <- 1000
iters <- 1000
data_length <- c(50,100,200,500,1000)
gamma0 <- 0.5

num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)

clusterEvalQ(cl, {
  library(Hmisc)
  library(nnet)
})

clusterExport(cl, c(
  "drm_1", "emp", "CEL_est", "ReadBasis", "read.DATA", "quan.ecdf",
  "quantiles", "BasisFun", "B", "gamma0", "wtd.quantile"
))

for(n in data_length){
  cat("Processing n =", n, "on", num_cores, "cores...\n")
  
  ### -------------------------
  ### dim1: gamma_real=0.5
  sim_dim1 <- function(iter, n) {
    set.seed(-iter)
    x1 <- 0.5*rnorm(n, mean=19, sd=0.5)+0.5*rnorm(n, mean=21, sd=0.5)
    y1 <- 0.5*rnorm(n, mean=17, sd=2)+0.5*rnorm(n, mean=23, sd=2)
    DATA <- list(cbind(1:n, x1), cbind(1:n, y1))
    gamma_hat_drm <- drm_1(DATA, quantiles, BasisFun, c(1,2))
    gamma_hat_emp <- emp(DATA, c(1,2))
    
    gamma_hat_drm_star <- numeric(B)
    gamma_hat_emp_star <- numeric(B)
    for (k in 1:B) {
      set.seed(B * iter + k)
      boot_DATA <- list(
        cbind(1:n, sample(x1, n, replace = TRUE)),
        cbind(1:n, sample(y1, n, replace = TRUE))
      )
      gamma_hat_drm_star[k] <- drm_1(boot_DATA, quantiles, BasisFun, c(1,2))
      gamma_hat_emp_star[k] <- emp(boot_DATA, c(1,2))
    }
    
    temp_drm <- mean(gamma_hat_drm_star - gamma_hat_drm > gamma_hat_drm - gamma0)
    p_value_drm <- 2 * min(temp_drm, 1 - temp_drm)
    temp_emp <- mean(gamma_hat_emp_star - gamma_hat_emp > gamma_hat_emp - gamma0)
    p_value_emp <- 2 * min(temp_emp, 1 - temp_emp)
    return(list(p_value_drm = p_value_drm, p_value_emp = p_value_emp))
  }
  
  result <- parSapply(cl, 1:iters, sim_dim1, n = n)
  p_value_drm_dim1 <- as.numeric(result["p_value_drm",])
  p_value_emp_dim1 <- as.numeric(result["p_value_emp",])
  cat("Finished dim1 for n =", n, "\n")
  
  ### -------------------------
  ### dim2: gamma_real=0.3097
  sim_dim2 <- function(iter, n) {
    set.seed(iter)
    x2 <- 0.5*rnorm(n, mean=19, sd=0.5)+0.5*rnorm(n, mean=21, sd=0.5)
    y2 <- 0.34*rnorm(n, mean=17, sd=2)+0.66*rnorm(n, mean=23, sd=2)
    DATA <- list(cbind(1:n, x2), cbind(1:n, y2))
    gamma_hat_drm <- drm_1(DATA, quantiles, BasisFun, c(1,2))
    gamma_hat_emp <- emp(DATA, c(1,2))
    
    gamma_hat_drm_star <- numeric(B)
    gamma_hat_emp_star <- numeric(B)
    for (k in 1:B) {
      set.seed(-B * iter - k)
      boot_DATA <- list(
        cbind(1:n, sample(x2, n, replace = TRUE)),
        cbind(1:n, sample(y2, n, replace = TRUE))
      )
      gamma_hat_drm_star[k] <- drm_1(boot_DATA, quantiles, BasisFun, c(1,2))
      gamma_hat_emp_star[k] <- emp(boot_DATA, c(1,2))
    }
    
    temp_drm <- mean(gamma_hat_drm_star - gamma_hat_drm > gamma_hat_drm - gamma0)
    p_value_drm <- 2 * min(temp_drm, 1 - temp_drm)
    temp_emp <- mean(gamma_hat_emp_star - gamma_hat_emp > gamma_hat_emp - gamma0)
    p_value_emp <- 2 * min(temp_emp, 1 - temp_emp)
    return(list(p_value_drm = p_value_drm, p_value_emp = p_value_emp))
  }
  
  result <- parSapply(cl, 1:iters, sim_dim2, n = n)
  p_value_drm_dim2 <- as.numeric(result["p_value_drm",])
  p_value_emp_dim2 <- as.numeric(result["p_value_emp",])
  cat("Finished dim2 for n =", n, "\n")
  
  p_value_drm <- c(p_value_drm_dim1, p_value_drm_dim2)
  adjusted_p_drm <- p.adjust(p_value_drm, method = "BH")
  p_value_emp <- c(p_value_emp_dim1, p_value_emp_dim2)
  adjusted_p_emp <- p.adjust(p_value_emp, method = "BH")
  
  write.csv(p_value_drm, file = paste0("p_value_drm_twocompnormal_", n, ".csv"), row.names = FALSE)
  write.csv(adjusted_p_drm, file = paste0("adjusted_p_value_drm_twocompnormal_", n, ".csv"), row.names = FALSE)
  write.csv(p_value_emp, file = paste0("p_value_emp_twocompnormal_", n, ".csv"), row.names = FALSE)
  write.csv(adjusted_p_emp, file = paste0("adjusted_p_value_emp_twocompnormal_", n, ".csv"), row.names = FALSE)
}

stopCluster(cl)
cat("All simulations finished!\n")

