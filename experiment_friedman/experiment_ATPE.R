#### Friedman Test + Nemenyi Post-hoc (Demšar, 2006)
library(readxl)

### 1. Input
s <- 7  # number of algorithms
n_datasets <- 50  # number of datasets

# Algorithm names (for readability)
alg_names <- c("SegRNN_96_96_ATPE", "ModernTCN_96_96_ATPE", "iTransformer_96_96_ATPE",
               "MTST_96_96_ATPE", "DLinear_96_96_ATPE", "NLinear_96_96_ATPE", "MTSMixer_96_96_ATPE")

cols_n_rows <- list(Row = alg_names, Col = alg_names)

### 2. Data Preparation
data <- read_excel("results_all.xlsx")

# Extract each algorithm's ATPE vector (lower is better)
x0 <- as.vector(as.matrix(data['SegRNN_96_96_ATPE']))
x1 <- as.vector(as.matrix(data['ModernTCN_96_96_ATPE']))
x2 <- as.vector(as.matrix(data['iTransformer_96_96_ATPE']))
x3 <- as.vector(as.matrix(data['MTST_96_96_ATPE']))
x4 <- as.vector(as.matrix(data['DLinear_96_96_ATPE']))
x5 <- as.vector(as.matrix(data['NLinear_96_96_ATPE']))
x6 <- as.vector(as.matrix(data['MTSMixer_96_96_ATPE']))

X_list <- list(x0, x1, x2, x3, x4, x5, x6)

### 3. Step 1: Compute ranks per dataset
# Create a matrix: rows = datasets, cols = algorithms
perf_matrix <- do.call(cbind, X_list)  # n_datasets x s

# Since lower ATPE is better, rank from low to high (best = rank 1)
ranks_matrix <- t(apply(perf_matrix, 1, rank, ties.method = "average"))

# Step 2: Compute average ranks for each algorithm
avg_ranks <- colMeans(ranks_matrix)

### 4. Step 3: Friedman test (overall significance)
# We'll use the F-distribution approximation (Iman & Davenport)
chi2_F <- (12 * n_datasets) / (s * (s + 1)) * (sum(avg_ranks^2) - s * (s + 1)^2 / 4)
F_F <- (n_datasets - 1) * chi2_F / (n_datasets * (s - 1) - chi2_F)
p_friedman <- pf(F_F, df1 = s - 1, df2 = (s - 1) * (n_datasets - 1), lower.tail = FALSE)

cat("Friedman test p-value:", p_friedman, "\n")

### 5. Step 4: Nemenyi post-hoc test
# Critical value q_alpha for alpha = 0.05, k = s algorithms
# Use studentized range distribution (qtukey)
alpha <- 0.05
q_alpha <- qtukey(1 - alpha, nmeans = s, df = Inf)  # df=Inf approximates large N

# Critical Difference (CD)
CD <- q_alpha / sqrt(2) * sqrt(s * (s + 1) / (6 * n_datasets))

cat("Critical Difference (CD):", round(CD, 4), "\n")

### 6. Initialize output matrices (lower triangle only)
rank_diff <- array(NA, dim = c(s, s))
significance <- array(0, dim = c(s, s))  # 0 = not significant

dimnames(rank_diff) <- cols_n_rows
dimnames(significance) <- cols_n_rows

# Fill lower triangle (i > j)
for (i in 2:s) {
  for (j in 1:(i - 1)) {
    diff_ij <- avg_ranks[i] - avg_ranks[j]
    rank_diff[i, j] <- diff_ij
    
    if (abs(diff_ij) > CD) {
      # If row algorithm has higher rank -> worse performance
      if (diff_ij > 0) {
        significance[i, j] <- -1  # row algo significantly WORSE than column
      } else {
        significance[i, j] <- 1   # row algo significantly BETTER than column
      }
    } else {
      significance[i, j] <- 0     # not significant
    }
  }
}

### 7. Save results

# (1) Average ranks (vector)
write.csv(data.frame(AvgRank = avg_ranks), file = "ATPE_friedman_avg_rank.csv", row.names = TRUE)

# (2) Rank differences (lower triangle) + CD in filename
CD_str <- ifelse(CD < 0.0001, "0.0000", sprintf("%.4f", round(CD, 4)))
write.csv(rank_diff, file = paste0("ATPE_friedman_rank_diff_CD_", CD_str, ".csv"))

# (3) Significance direction (lower triangle)
write.csv(significance, file = paste0("ATPE_friedman_significance_direction_CD_", CD_str, ".csv"))

cat("Friedman + Nemenyi results saved.\n")