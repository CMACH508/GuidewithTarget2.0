#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

# Check for correct number of arguments
if (length(args) < 9) {
  stop("Usage: Rscript script_name.R input.txt output.bed cleavage_group non_cleavage_group num_cleavage_reps num_non_cleavage_reps [logfc_thres|nofilter] pval_thres Y")
}

input_file <- args[1]
output_file <- args[2]
cleavage_group <- args[3]
non_cleavage_group <- args[4]
num_cleavage_reps <- as.integer(args[5])
num_non_cleavage_reps <- as.integer(args[6])

# --- 修改点 1: 处理 nofilter 逻辑 ---
filter_flag <- TRUE
if (args[7] == "nofilter") {
  filter_flag <- FALSE
  logfc_thres <- NA
  pval_thres <- NA
  cat("!!! Detected 'nofilter' mode: All rows will be retained regardless of thresholds !!!\n")
} else {
  logfc_thres <- as.numeric(args[7])
  pval_thres <- as.numeric(args[8])
}
keep_zeroinMUT <- args[9]

# Print input parameters for user confirmation
cat("=== INPUT PARAMETERS ===\n")
cat(paste("Input file:", input_file, "\n"))
cat(paste("Output file:", output_file, "\n"))
cat(paste("Log2FC threshold:", args[7], "\n"))
cat(paste("P-value threshold:", args[8], "\n"))
cat(paste("if discard non-zero rows in MUT:", keep_zeroinMUT, "\n"))
cat("========================\n\n")

# Read the data
data <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Generate column names
cleavage_cols <- paste0(cleavage_group, "_", 1:num_cleavage_reps, "_norm")
non_cleavage_cols <- paste0(non_cleavage_group, "_", 1:num_non_cleavage_reps, "_norm")

# Check if all required columns exist
required_cols <- c(cleavage_cols, non_cleavage_cols)
missing_cols <- setdiff(required_cols, colnames(data))
if (length(missing_cols) > 0) {
  stop(paste("Missing columns in input file:", paste(missing_cols, collapse = ", ")))
}

# Perform statistical analysis
if(length(cleavage_cols) > 1){
  mean_cleavage <- rowMeans(data[, cleavage_cols], na.rm = TRUE)
  mean_non_cleavage <- rowMeans(data[, non_cleavage_cols], na.rm = TRUE)
}else{
  mean_cleavage <- data[, cleavage_cols]
  mean_non_cleavage <- data[, non_cleavage_cols]
}

# Compute log2 fold change
data$log2FC <- log2(mean_non_cleavage / mean_cleavage)

# --- 修改点 2: 分支处理过滤逻辑 ---
if (num_cleavage_reps == 1 && num_non_cleavage_reps == 1) {
  cat("Single replicate case detected.\n")
  data$p_value <- "."
  
  if (filter_flag) {
    sig_down <- data[which(data$log2FC < logfc_thres & !is.na(data$log2FC)), ]
  } else {
    sig_down <- data
  }
  
} else {
  # Multiple replicates case
  has_na <- apply(data[, required_cols], 1, function(x) any(is.na(x)))
  ss_cleavage <- rowSums((data[, cleavage_cols] - mean_cleavage)^2)
  ss_non_cleavage <- rowSums((data[, non_cleavage_cols] - mean_non_cleavage)^2)
  
  pooled_var <- (ss_cleavage + ss_non_cleavage) / (num_cleavage_reps + num_non_cleavage_reps - 2)
  se <- sqrt(pooled_var * (1/num_cleavage_reps + 1/num_non_cleavage_reps))
  
  t_stat <- (mean_non_cleavage - mean_cleavage) / se
  t_stat[is.nan(t_stat)] <- 0
  
  df <- num_cleavage_reps + num_non_cleavage_reps - 2
  data$p_value <- 2 * pt(-abs(t_stat), df = df)
  data$p_value[has_na] <- NA
  data$log2FC[has_na] <- NA
  
  if (filter_flag) {
    sig_down <- data[which(data$p_value < pval_thres & data$log2FC < logfc_thres & !is.na(data$p_value)), ]
  } else {
    sig_down <- data
  }
}

# --- 修改点 3: 只有在 filter_flag 为 TRUE 时才执行 keep_zeroinMUT 过滤 ---
if (filter_flag && keep_zeroinMUT == "Y") {
  cat("Filtering: Keeping only rows where all", non_cleavage_group, "replicates are zero\n")
  # 注意：这里需要根据原始数据的列名规则匹配
  curr_non_cleavage_cols <- grep(paste0("^", non_cleavage_group, "_[0-9]+_norm$"), colnames(sig_down), value = TRUE)
  sig_down <- sig_down[rowSums(sig_down[, curr_non_cleavage_cols, drop = FALSE] == 0, na.rm = TRUE) == length(curr_non_cleavage_cols), ]
  cat("After zero-filtering,", nrow(sig_down), "rows remain\n")
} else if (!filter_flag) {
  cat("No-filter mode: Skipping zero-value filtering.\n")
}

# Create BED file content
# 增加处理：如果是 nofilter 模式，可能会有 NA 值，需要剔除无法转为 BED 的行（如没有 chr 的行）
sig_down <- sig_down[!is.na(sig_down$chr), ]

bed_data <- data.frame(
  chr = sig_down$chr,
  start = sig_down$start,
  end = sig_down$end,
  name = sig_down$uniqpos,
  score = "all",
  strand = sig_down$strand
)
bed_data <- unique(bed_data)

# Write to BED file
write.table(bed_data, file = output_file, sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

cat(paste("\nAnalysis complete. Total rows in output:", nrow(bed_data), "\n"))