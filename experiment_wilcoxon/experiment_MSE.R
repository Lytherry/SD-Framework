#### pairwise Wilcoxon signed-rank tests with the BH correction
library(readxl)

### 1. Input
s <- 7 # number of algorithms
cols_n_rows <- list(
  Row = c("SegRNN_96_96_mse","ModernTCN_96_96_mse","iTransformer_96_96_mse",
          "MTST_96_96_mse","DLinear_96_96_mse","NLinear_96_96_mse","MTSMixer_96_96_mse"), 
  Col = c("SegRNN_96_96_mse","ModernTCN_96_96_mse","iTransformer_96_96_mse",
          "MTST_96_96_mse","DLinear_96_96_mse","NLinear_96_96_mse","MTSMixer_96_96_mse")
)

### 2. Data Preparation
data <- read_excel("results_all.xlsx")

x0 <- as.vector(as.matrix(data['SegRNN_96_96_mse']))
x1 <- as.vector(as.matrix(data['ModernTCN_96_96_mse']))
x2 <- as.vector(as.matrix(data['iTransformer_96_96_mse']))
x3 <- as.vector(as.matrix(data['MTST_96_96_mse']))
x4 <- as.vector(as.matrix(data['DLinear_96_96_mse']))
x5 <- as.vector(as.matrix(data['NLinear_96_96_mse']))
x6 <- as.vector(as.matrix(data['MTSMixer_96_96_mse']))

X_list <- list(x0, x1, x2, x3, x4, x5, x6)  # store in a list for easier looping

### 3. Initialize output arrays with NA
median_diff <- array(NA, dim = c(s, s))
p_value_wilcox <- array(NA, dim = c(s, s))

dimnames(median_diff) <- cols_n_rows
dimnames(p_value_wilcox) <- cols_n_rows

### 4. Pairwise Wilcoxon signed-rank tests
for (i in 2:s) {
  for (j in 1:(i - 1)) {
    # Perform two-sided Wilcoxon signed-rank test
    test_res <- wilcox.test(X_list[[i]], X_list[[j]], 
                            paired = TRUE, 
                            alternative = "two.sided",
                            exact = FALSE,
                            correct = TRUE)
    
    # Store median difference: algo_i - algo_j
    median_diff[i, j] <- median(X_list[[i]] - X_list[[j]])
    
    # Store p-value
    p_value_wilcox[i, j] <- test_res$p.value
  }
}

### 5. Save results
write.csv(median_diff, file = "mse_wilcox_median_diff.csv")
write.csv(p_value_wilcox, file = "mse_wilcox_p_value.csv")

cat("Wilcoxon signed-rank test results saved.\n")