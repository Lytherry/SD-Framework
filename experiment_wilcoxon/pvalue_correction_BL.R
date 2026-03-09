# ----------------------
# 步骤1：读取所有CSV文件并保留行列名
# ----------------------
input_files <- c("mse_wilcox_p_value.csv", "ATPE_wilcox_p_value.csv", "AIT_wilcox_p_value.csv")
output_files <- c("adjusted_mse_wilcox_p_value_BL.csv", "adjusted_ATPE_wilcox_p_value_BL.csv", "adjusted_AIT_wilcox_p_value_BL.csv")

# 读取三个文件为列表，保留行列名
matrix_list <- lapply(input_files, function(f) {
  read.csv(f, header = TRUE, row.names = 1)
})

# ----------------------
# 步骤2：合并所有p值并记录位置信息
# ----------------------
# 创建数据框记录每个p值的位置（文件ID, 行名, 列名, 原始位置索引）
p_value_info <- data.frame()

for (file_id in seq_along(matrix_list)) {
  mat <- matrix_list[[file_id]]
  rownames <- rownames(mat)
  colnames <- colnames(mat)
  
  # 遍历矩阵每个位置
  for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
      p <- mat[i, j]
      if (!is.na(p)) {
        # 记录非NA值的位置信息
        p_value_info <- rbind(p_value_info, data.frame(
          file_id = file_id,
          row = rownames[i],
          col = colnames[j],
          row_index = i,
          col_index = j,
          p_value = p
        ))
      }
    }
  }
}

# ----------------------
# 步骤3：整体进行BL校正
# ----------------------
p_adjust_BL <- function(p_values) {
  n <- length(p_values)
  if (n <= 0) stop("p_values不能为空")
  
  # 计算调和数 H(n) = 1 + 1/2 + ... + 1/n
  Hn <- sum(1 / seq_len(n))
  
  # 排序处理
  o <- order(p_values)
  sorted_p <- p_values[o]
  
  # BL调整公式
  adjusted <- sorted_p * n * Hn / seq_len(n)
  
  # 保证单调非递增
  adjusted <- cummin(adjusted)
  
  # 限制在[0,1]区间
  adjusted <- pmin(adjusted, 1)
  
  # 恢复原始顺序
  adjusted[order(o)]
}

all_p <- p_value_info$p_value
adjusted_p <- p_adjust_BL(all_p)
p_value_info$adjusted_p <- adjusted_p

# ----------------------
# 步骤4：将校正值填回原始矩阵结构
# ----------------------
# 创建新的校正矩阵列表
adjusted_list <- lapply(matrix_list, function(mat) {
  # 初始化与原始矩阵相同结构的NA矩阵
  corrected_mat <- matrix(NA, 
                          nrow = nrow(mat),
                          ncol = ncol(mat),
                          dimnames = list(rownames(mat), colnames(mat)))
  corrected_mat
})

# 填充校正后的值
for (k in 1:nrow(p_value_info)) {
  info <- p_value_info[k, ]
  # 通过行列索引定位到具体矩阵位置
  adjusted_list[[info$file_id]][info$row_index, info$col_index] <- info$adjusted_p
}

# ----------------------
# 步骤5：保存带行列名的CSV文件
# ----------------------
for (file_id in seq_along(adjusted_list)) {
  write.csv(adjusted_list[[file_id]], output_files[file_id])
}