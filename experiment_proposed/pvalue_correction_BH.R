# ----------------------
# step 1：read .csv files
# ----------------------
input_files <- c("mse_p_value.csv", "ATPE_p_value.csv", "AIT_p_value.csv")
output_files <- c("adjusted_mse_p_value_BH.csv", "adjusted_ATPE_p_value_BH.csv", "adjusted_AIT_p_value_BH.csv")

matrix_list <- lapply(input_files, function(f) {
  read.csv(f, header = TRUE, row.names = 1)
})

# ----------------------
# step 2：merge all p-values and record the location information
# ----------------------
p_value_info <- data.frame()

for (file_id in seq_along(matrix_list)) {
  mat <- matrix_list[[file_id]]
  rownames <- rownames(mat)
  colnames <- colnames(mat)
  
  for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
      p <- mat[i, j]
      if (!is.na(p)) {
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
# step 3：Perform BH correction on the entire dataset
# ----------------------
all_p <- p_value_info$p_value
adjusted_p <- p.adjust(all_p, method = "BH")

p_value_info$adjusted_p <- adjusted_p

# ----------------------
# step 4：fill the corrected values back into the original matrix structure
# ----------------------
adjusted_list <- lapply(matrix_list, function(mat) {
  corrected_mat <- matrix(NA, 
                          nrow = nrow(mat),
                          ncol = ncol(mat),
                          dimnames = list(rownames(mat), colnames(mat)))
  corrected_mat
})

for (k in 1:nrow(p_value_info)) {
  info <- p_value_info[k, ]
  adjusted_list[[info$file_id]][info$row_index, info$col_index] <- info$adjusted_p
}

# ----------------------
# step 5：save files
# ----------------------
for (file_id in seq_along(adjusted_list)) {
  write.csv(adjusted_list[[file_id]], output_files[file_id])
}