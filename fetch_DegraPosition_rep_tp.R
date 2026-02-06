args <- commandArgs(trailingOnly = TRUE)
cleav_suffix <- args[1] # WT
nocleav_suffix <- args[2] # KO
data_add <- args[3]
num_of_rep_cleav <- as.numeric(args[4]) # 3
num_of_rep_nocleav <- as.numeric(args[5])
set_range <- as.numeric(args[6]) # 50
risc <- args[7]
targetseq_len <- args[8]

# cleav_suffix = "RaDeWT"
# nocleav_suffix = "RaDeMut"
# data_add = "/data/zenglin/rabbit/degradome/wkdir/PIWIL3_PIWIL3_associated_RaDeWT_RaDeMU/"
# num_of_rep_cleav = 1
# num_of_rep_nocleav = 1
# set_range = 20
# risc = "PIWIL3_PIWIL3_associated"
# targetseq_len = "/data/zenglin/rabbit/degradome/transcriptome/OryCun_seq_tp_stat.txt"


seqlength <- read.table(targetseq_len, sep = "\t", header = FALSE)
colnames(seqlength) <- c("seqid", "seqlen")

options(scipen = 100)
library("data.table")
suppressPackageStartupMessages(library("dplyr"))
total_sample <- num_of_rep_cleav + num_of_rep_nocleav

#risc = "ago2_guide1"

#flagstat1 <- readLines(paste0(cleav_suffix,"_R2.trimmed.flagstat.txt"))
#flagstat_cleav <- as.numeric(unlist(strsplit(flagstat1[5]," "))[1])

#flagstat2 <- readLines(paste0(nocleav_suffix,"_R2.trimmed.flagstat.txt"))
#flagstat_nocleav <- as.numeric(unlist(strsplit(flagstat2[5]," "))[1]) # nolint


# Define the function
find_matched_index <- function(flagstat) {
  pattern <- "^\\d+\\s*\\+\\s*0\\s*mapped\\s*\\(\\d+\\.\\d+%\\s*:\\s*N/A\\)$"
  matched_indices <- grep(pattern, flagstat)
  valid_indices <- matched_indices[sapply(matched_indices, function(idx) {
    # Split the string to extract the first integer part
    parts <- strsplit(flagstat[idx], " ")[[1]]
    as.integer(parts[1]) != 0
  })]
  return(valid_indices)
}

process_cleav <- function(cleav_suffix_rep) {
  print(paste0("Processing cleavage: ", cleav_suffix_rep))
  flagstat1 <- readLines(paste0(cleav_suffix_rep, "_transposon.flagstat.txt"))
  matched_index <- find_matched_index(flagstat1)
  flagstat_cleav <- as.numeric(
    unlist(strsplit(flagstat1[matched_index], " "))[1]
  )  # [7]
  cleav_samfile <- paste0(
    data_add,
    "/",
    cleav_suffix_rep,
    "_transposon.candidates.sam"
  )
  cleav <- fread(cleav_samfile, sep = "\t")
  col2rename <- c("id", "strand", "chr", "pos_1_based", "sequence")
  colnames(cleav) <- col2rename
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!bed file 0-based  not include last character
  cleav$strand <- ifelse(cleav$strand == 0, "+", "-")
  range2ext <- set_range
  #if(reference=="1"){
  #	cleav <- cleav[which(cleav$strand=="+"),]
  #	range2ext = set_range
  #}
  cleav$start <- as.numeric(cleav$pos_1_based) - range2ext - 1
  cleav$start <- ifelse(cleav$start < 0, 0, cleav$start)
  cleav$chr <- gsub("#", "_", cleav$chr)
  cleav <- merge(cleav, seqlength, by.x = "chr", by.y = "seqid")
  cleav$tmp <- as.numeric(cleav$pos_1_based) + range2ext
  cleav$end <- pmin(cleav$tmp, cleav$seqlen)
  cleav$tmp <- NULL
  cleav$seqlen <- NULL
  # cleav$id <- unlist(strsplit(cleav$id,"\\:"))[seq(2,nrow(cleav)*2,2)]
  cleav$uniqpos <- paste0(cleav$chr, "_", cleav$pos_1_based, "_", cleav$strand)
  products_stat <- as.data.frame(table(cleav$uniqpos))
  cleav <- merge(cleav, products_stat, by.x = "uniqpos", by.y = "Var1")
  cleav <- cleav %>%
    select(uniqpos, id, strand, chr, pos_1_based, start, end, Freq) %>%
    group_by(uniqpos) %>%
    mutate(combi_id = paste(id, collapse = " | "))
  cleav <- unique(as.data.frame(
    cleav[, c("chr", "start", "end", "uniqpos", "Freq",
              "strand", "combi_id", "pos_1_based")]
  ))
  write.table(
    file = paste0(
      data_add,
      cleav_suffix_rep,
      "_3CleavProcutsStatRaw.tp.txt"
    ),
    cleav,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  reduced_cleav <- cleav[, c("chr", "start", "end",
                             "uniqpos", "Freq", "strand", "pos_1_based")]
  reduced_cleav$NormFreq <- round(
    (as.numeric(reduced_cleav$Freq) + 1) / flagstat_cleav * 1e6, 4
  )
  psuedo_cnt <- round((1 / flagstat_cleav * 1e6), 4)
  reduced_cleav <- reduced_cleav[order(reduced_cleav$Freq, decreasing = TRUE), ]
  colnames(reduced_cleav) <- c(
    "chr", "start", "end", "uniqpos", cleav_suffix_rep, "strand",
    "pos_1_based", paste0(cleav_suffix_rep, "_norm")
  )
  write.table(
    file = paste0(
      data_add,
      cleav_suffix_rep,
      "_3CleavProductsStatNorm.tp.txt"
    ),
    reduced_cleav,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  reduced_cleav$source <- cleav_suffix_rep
  return(list(df = reduced_cleav, cnt = psuedo_cnt))
}



# Initiate container to store all contents finally
cnt_vct <- c()
events <- data.frame()


# For Cleavage
if (num_of_rep_cleav == 1) {
  rep_name <- paste0(cleav_suffix, "_1")
  result <- process_cleav(rep_name)
  tmp <- result$df[, c(4, 5, 8)]
  cnt_vct <- c(cnt_vct, result$cnt)
  colnames(result$df)[5] <- "events_freq"
  colnames(result$df)[8] <- "events_norm"
  events <- rbind(events, result$df)
}else {
  rep_name <- paste0(cleav_suffix, "_1") # _1 or 1
  result <- process_cleav(rep_name)
  tmp <- result$df[, c(4, 5, 8)]
  cnt <- result$cnt
  cnt_vct <- c(cnt_vct, cnt)
  colnames(result$df)[5] <- "events_freq"
  colnames(result$df)[8] <- "events_norm"
  events <- rbind(events, result$df)
  for (nrep in c(2:num_of_rep_cleav)){
    rep_name <- paste0(cleav_suffix, "_", nrep)
    result <- process_cleav(rep_name)
    tmp0 <- result$df
    cnt0 <- result$cnt
    tmp <- merge(tmp, tmp0[, c(4, 5, 8)], by = "uniqpos", all = TRUE)
    cnt_vct <- c(cnt_vct, cnt0)
    colnames(result$df)[5] <- "events_freq"
    colnames(result$df)[8] <- "events_norm"
    events <- rbind(events, result$df)
  }
  #for(col_idx in 2:ncol(tmp)){
  #	# tmp <- tmp %>% na.fill(cnt_vct[col_idx], .cols = 1)	
  #	tmp[is.na(tmp[, col_idx]), col_idx] <- cnt_vct[col_idx-1]
  #}
}
cleav_tmp <- tmp


# For Nocleavage
if (num_of_rep_nocleav == 1) {
  rep_name <- paste0(nocleav_suffix, "_1")
  result <- process_cleav(rep_name)
  tmp <- result$df[, c(4, 5, 8)]
  cnt_vct <- c(cnt_vct, result$cnt)
  colnames(result$df)[5] <- "events_freq"
  colnames(result$df)[8] <- "events_norm"
  events <- rbind(events, result$df)
}else {
  rep_name <- paste0(nocleav_suffix, "_1")
  result <- process_cleav(rep_name)
  tmp <- result$df[, c(4, 5, 8)]
  cnt <- result$cnt
  cnt_vct <- c(cnt_vct, cnt)
  colnames(result$df)[5] <- "events_freq"
  colnames(result$df)[8] <- "events_norm"
  events <- rbind(events, result$df)
  for (nrep in c(2:num_of_rep_nocleav)){
    rep_name <- paste0(nocleav_suffix, "_", nrep)
    result <- process_cleav(rep_name)
    tmp0 <- result$df
    cnt0 <- result$cnt
    tmp <- merge(tmp, tmp0[, c(4, 5, 8)], by = "uniqpos", all = TRUE)
    cnt_vct <- c(cnt_vct, cnt0)
    colnames(result$df)[5] = "events_freq"
    colnames(result$df)[8] = "events_norm"
    events <- rbind(events, result$df)
  }
  #for(col_idx in 2:ncol(tmp)){
  #       # tmp <- tmp %>% na.fill(cnt_vct[col_idx], .cols = 1)
  #       tmp[is.na(tmp[, col_idx]), col_idx] <- cnt_vct[col_idx-1]
  #}

}
nocleav_tmp <- tmp

total_table <- merge(cleav_tmp, nocleav_tmp, by = "uniqpos", all = TRUE)
reorder <- c(paste0(cleav_suffix, "_", c(1:num_of_rep_cleav)),
             paste0(nocleav_suffix, "_",c(1:num_of_rep_nocleav)),
             paste0(cleav_suffix, "_", c(1:num_of_rep_cleav), "_norm"),
             paste0(nocleav_suffix, "_", c(1:num_of_rep_nocleav), "_norm"))
#total_table[,2:ncol(total_table)] <- total_table[,reorder]
total_table <- total_table[, c("uniqpos", reorder)]

# FIXME
# Separate raw counts and normalized counts clearly
raw_count_cols <- 2:(total_sample + 1)
norm_count_cols <- (total_sample + 2):ncol(total_table)

# For RAW COUNTS: use 0 for missing values (true absence)
total_table[, raw_count_cols] <- lapply(total_table[, raw_count_cols], 
                                       function(x) replace(x, is.na(x), 0))
# For NORMALIZED COUNTS: use appropriate pseudo-counts
for (i in seq_along(norm_count_cols)) {
  col_idx <- norm_count_cols[i]
  total_table[is.na(total_table[, col_idx]), col_idx] <- cnt_vct[i]
}

idx <- 1
for (col_idx in (total_sample + 2):ncol(total_table)){
       print(cnt_vct[idx])
       total_table[is.na(total_table[, col_idx]), col_idx] <- cnt_vct[idx]
       idx <- idx + 1
}



if (FALSE) {
  if (total_sample == 2) {
    start_group1 <- total_sample + 2
    end_group1 <- total_sample + 2 + num_of_rep_cleav - 1
    start_group2 <- num_of_rep * 2 + 2 + num_of_rep
    end_group2 <- num_of_rep * 2 + 2 + 2 * num_of_rep - 1
    total_table$mean_group1 <- c(total_table[, start_group1:end_group1])
    total_table$mean_group2 <- c(total_table[, start_group2:end_group2])
    total_table$mean_ratio <- total_table$mean_group1 / total_table$mean_group2
  }else {
    start1 <- num_of_rep * 2 + 2
    end1 <- start1 + num_of_rep - 1
    start2 <- end1 + 1
    end2 <- start2 + num_of_rep - 1
    total_table$mean_group1 <- rowMeans(
      total_table[, start1:end1],
      na.rm = TRUE
    )
    total_table$mean_group2 <- rowMeans(
      total_table[, start2:end2],
      na.rm = TRUE
    )
    total_table$mean_ratio <- total_table$mean_group1 / total_table$mean_group2
  }
}

# Perform t-test row-wise and calculate p-values
#p_values <- apply(total_table[,((num_of_rep*2+2)):ncol(total_table)], 1, function(x) {
#  group1 <- x[1:num_of_rep]
#  group2 <- x[(num_of_rep+1):(2*(num_of_rep))]
#  t_test_result <- t.test(group1, group2, paired = TRUE)
#  return(t_test_result$p.value)
#})

# Add p-values to the dataframe
#total_table$p_value <- p_values
#sig_total_table <- total_table[which(total_table$p_value<0.05 & total_table$mean_ratio>1),]



events$placeholder <- events$source
events$start <- as.numeric(events$start)
events$end <- as.numeric(events$end)
columns_to_select <- c("chr", "start", "end", "uniqpos", "strand",
                       "pos_1_based", "placeholder")
events <- unique(events[, columns_to_select])
events <- merge(events, total_table, by = "uniqpos")
events_uniq <- events
events_uniq$placeholder <- "all"
events_uniq <- unique(events_uniq)
events_uniq$placeholder <- paste0("products_",1:nrow(events_uniq))


output_file <- paste0(data_add, risc, "_insert_expr.tp.txt")
write.table(
  file = output_file,
  events_uniq,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

output_file <- paste0(data_add, risc, "_insert.tp.bed")
selected_columns <- c("chr", "start", "end", "uniqpos", "placeholder", "strand")
write.table(
  file = output_file,
  events_uniq[, selected_columns],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)