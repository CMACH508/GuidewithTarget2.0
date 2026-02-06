args = commandArgs(T)
file1 = args[1] # POI_Target_Pair:ago3_POI_Target_Pair.txt
risc = as.character(args[2])
options(scipen=999)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library(parallel))

# Enhanced logging system - define before any parallel operations
log_file <- paste0(risc, "_debug.log")
debug_log <- function(msg, level = "INFO") {
  timestamp_msg <- paste(Sys.time(), level, msg)
  write(timestamp_msg, file = log_file, append = TRUE)
  print(timestamp_msg)
}

debug_log(paste("Starting processing of", file1))

# Safe file reading
result <- tryCatch({
  lines <- readLines(paste0(file1))
  debug_log(paste("Successfully read file with", length(lines), "lines"))
  lines
}, error = function(e) {
  debug_log(paste("FATAL: Error reading file", e$message), "ERROR")
  stop("File reading failed")
})

rownum_query <- which(startsWith(result,"Query="))
debug_log(paste("Found", length(rownum_query), "queries"))

if (length(rownum_query) == 0) {
  debug_log("No queries found in file!", "ERROR")
  stop("No queries found")
}

# Create a simplified debug function for parallel workers
worker_debug_log <- function(msg, level = "INFO") {
  timestamp_msg <- paste(Sys.time(), Sys.getpid(), level, msg)
  # Workers write to a separate file to avoid conflicts
  worker_log_file <- paste0(risc, "_worker_debug.log")
  write(timestamp_msg, file = worker_log_file, append = TRUE)
}

# Enhanced ParseEachQuery with detailed debugging
ParseEachQuery <- function(query_num){
  worker_debug_log(paste("Starting query", query_num), "DEBUG")
  
  tryCatch({
    # Validate input parameters
    if (is.na(query_num)) {
      worker_debug_log(paste("Query number is NA"), "ERROR")
      return(NULL)
    }
    
    if (query_num < 1 || query_num > length(rownum_query)) {
      worker_debug_log(paste("Query number", query_num, "out of range [1,", length(rownum_query), "]"), "ERROR")
      return(NULL)
    }
    
    # Safe query extraction
    query_row <- rownum_query[query_num]
    worker_debug_log(paste("Query", query_num, "starts at row", query_row), "DEBUG")
    
    if (is.na(query_row) || query_row > length(result)) {
      worker_debug_log(paste("Invalid query_row:", query_row, "for query", query_num), "ERROR")
      return(NULL)
    }
    
    # Extract query block
    if(length(rownum_query) == 1){
      query_extrat <- result[query_row:length(result)]
    } else if(query_num == length(rownum_query)){
      query_extrat <- result[query_row:length(result)]
    } else {
      next_query <- rownum_query[query_num + 1]
      if (is.na(next_query) || next_query > length(result)) {
        query_extrat <- result[query_row:length(result)]
      } else {
        query_extrat <- result[query_row:(next_query - 1)]
      }
    }
    
    worker_debug_log(paste("Query", query_num, "block has", length(query_extrat), "lines"), "DEBUG")
    
    if (length(query_extrat) == 0) {
      worker_debug_log(paste("Empty query block for query", query_num), "ERROR")
      return(NULL)
    }
    
    # Initialize results dataframe
    df2 <- data.frame(query_id=character(), sub_id=character(), q_start=character(),
      q_end=character(), s_start=character(), s_end=character(),
      q_match_pos=character(), q_mismatch_pos=character(),
      identity=character(), alignment_length = numeric(), match_num = numeric(),
      mismatch_num = numeric(), stringsAsFactors=FALSE)
    
    # Extract basic info with debugging
    sub_id_list <- query_extrat[which(startsWith(query_extrat,">"))]
    worker_debug_log(paste("Query", query_num, "has", length(sub_id_list), "subjects"), "DEBUG")
    
    if (length(query_extrat) < 6) {
      worker_debug_log(paste("Query", query_num, "has insufficient lines:", length(query_extrat)), "ERROR")
      return(NULL)
    }
    
    query_id <- tryCatch({
      id <- substr(query_extrat[1], 8, nchar(query_extrat[1]))
      worker_debug_log(paste("Query ID:", id), "DEBUG")
      id
    }, error = function(e) {
      worker_debug_log(paste("Error extracting query ID:", e$message), "ERROR")
      paste0("query_", query_num)
    })
    
    Nohist <- query_extrat[6]
    worker_debug_log(paste("No hits status:", Nohist), "DEBUG")
    
    if(Nohist == "***** No hits found *****"){
      worker_debug_log(paste("No hits for query", query_num), "INFO")
    } else {
      # Process each subject
      for(subject in 1:length(sub_id_list)){
        subject_debug_prefix <- paste("Query", query_num, "Subject", subject)
        
        tryCatch({
          sub_id <- sub_id_list[subject]
          worker_debug_log(paste(subject_debug_prefix, "ID:", sub_id), "DEBUG")
          
          if (is.na(sub_id)) {
            worker_debug_log(paste(subject_debug_prefix, "has NA subject ID"), "ERROR")
            next
          }
          
          sub_rownum <- which(startsWith(query_extrat, sub_id))
          if (length(sub_rownum) == 0) {
            worker_debug_log(paste(subject_debug_prefix, "not found in query block"), "ERROR")
            next
          }
          
          sub_rownum1 <- min(sub_rownum + 11, length(query_extrat))
          worker_debug_log(paste(subject_debug_prefix, "rows:", sub_rownum, "to", sub_rownum1), "DEBUG")
          
          content_rel <- query_extrat[sub_rownum:sub_rownum1]
          
          if (length(content_rel) < 11) {
            worker_debug_log(paste(subject_debug_prefix, "insufficient content, length:", length(content_rel)), "ERROR")
            next
          }
          
          # Parse alignment details with error handling
          worker_debug_log(paste(subject_debug_prefix, "parsing alignment..."), "DEBUG")
          
          identity <- tryCatch({
            (unlist(strsplit(content_rel[5],"), Gaps = "))[1] %>% strsplit("\\(") %>% unlist)[2]
          }, error = function(e) {
            worker_debug_log(paste(subject_debug_prefix, "error parsing identity:", e$message), "ERROR")
            NA
          })
          
          alignment_summary <- tryCatch({
            (unlist(strsplit(content_rel[5]," Identities = "))[2] %>% strsplit(" \\(") %>% unlist)[1]
          }, error = function(e) {
            worker_debug_log(paste(subject_debug_prefix, "error parsing alignment summary:", e$message), "ERROR")
            NA
          })
          
          alignment_length <- tryCatch({
            as.numeric(unlist(strsplit(alignment_summary,"/"))[2])
          }, error = function(e) {
            worker_debug_log(paste(subject_debug_prefix, "error parsing alignment length:", e$message), "ERROR")
            NA
          })
          
          match_num <- tryCatch({
            as.numeric(unlist(strsplit(alignment_summary,"/"))[1])
          }, error = function(e) {
            worker_debug_log(paste(subject_debug_prefix, "error parsing match num:", e$message), "ERROR")
            NA
          })
          
          mismatch_num <- tryCatch({
            if (!is.na(alignment_length) && !is.na(match_num)) {
              alignment_length - match_num
            } else NA
          }, error = function(e) {
            worker_debug_log(paste(subject_debug_prefix, "error calculating mismatch num:", e$message), "ERROR")
            NA
          })
          
          # Parse positions
          temp1 <- tryCatch({
            temp <- unlist(strsplit(content_rel[8]," "))
            temp[which(temp != "")]
          }, error = function(e) {
            worker_debug_log(paste(subject_debug_prefix, "error splitting line 8:", e$message), "ERROR")
            character(0)
          })
          
          q_start <- if(length(temp1) >= 4) temp1[2] else NA
          q_end <- if(length(temp1) >= 4) temp1[4] else NA
          
          temp2 <- tryCatch({
            temp <- unlist(strsplit(content_rel[10]," "))
            temp[which(temp != "")]
          }, error = function(e) {
            worker_debug_log(paste(subject_debug_prefix, "error splitting line 10:", e$message), "ERROR")
            character(0)
          })
          
          s_end <- if(length(temp2) >= 4) temp2[2] else NA
          s_start <- if(length(temp2) >= 4) temp2[4] else NA
          
          worker_debug_log(paste(subject_debug_prefix, "positions - q_start:", q_start, "q_end:", q_end, "s_start:", s_start, "s_end:", s_end), "DEBUG")
          
          # Parse match positions
          matchpos <- content_rel[9]
          q_match_pos <- ""
          q_mismatch_pos <- ""
          
          if (!is.na(matchpos) && nchar(matchpos) > 0) {
            count_padding <- 0
            while(count_padding < nchar(matchpos) && substr(matchpos, count_padding+1, count_padding+1) != "|") {
              count_padding <- count_padding + 1
            }
            
            for(i in 1:nchar(matchpos)) {
              symbol <- substr(matchpos, i, i)
              if(symbol == "|") {
                q_match_pos <- c(q_match_pos, i - count_padding + as.numeric(q_start) - 1)
              }
            }
            q_match_pos <- q_match_pos[2:length(q_match_pos)]
            
            if (!is.na(q_start) && !is.na(q_end)) {
              q_mismatch_pos <- setdiff(c(as.numeric(q_start):(as.numeric(q_end) - 1)), q_match_pos)
            }
          }
          
          q_match_pos <- paste0(q_match_pos, collapse = ",")
          q_mismatch_pos <- paste0(q_mismatch_pos, collapse = ",")
          
          # Create result row
          pairresult <- data.frame(
            query_id = query_id, sub_id = sub_id, q_start = q_start, q_end = q_end,
            s_start = s_start, s_end = s_end, q_match_pos = q_match_pos,
            q_mismatch_pos = q_mismatch_pos, alignment_length = alignment_length,
            match_num = match_num, mismatch_num = mismatch_num, identity = identity
          )
          
          df2 <- rbind(df2, pairresult)
          worker_debug_log(paste(subject_debug_prefix, "processed successfully"), "DEBUG")
          
        }, error = function(e) {
          worker_debug_log(paste("ERROR in", subject_debug_prefix, ":", e$message), "ERROR")
        })
      }
    }
    
    # Post-process results
    if (nrow(df2) > 0) {
      worker_debug_log(paste("Query", query_num, "has", nrow(df2), "rows, starting post-processing"), "DEBUG")
      
      tryCatch({
        df2$dis1 <- 10 - as.numeric(df2$q_start)
        
        pattern_before_plus = ".*_(\\d+)_[+/-].*"
        df2$pos1_based <- as.numeric(gsub(pattern_before_plus, "\\1", df2$sub_id))
        
        pattern_before_parentheses = ".*-(\\d+)\\([+-]\\).*"
        df2$start_abs <- as.numeric(gsub(pattern_before_parentheses,"\\1", df2$sub_id))
        
        relpos <- abs(df2$pos1_based - df2$start_abs)
        dis <- relpos - as.numeric(df2$s_end)
        
        df2$rel_start <- as.numeric(df2$s_start) - relpos
        df2$rel_end <- as.numeric(df2$s_end) - relpos
        
        output_file <- paste0(risc, "_Parsed_NPOI_", query_num, ".txt")
        write.table(df2, file = output_file, row.names = FALSE, col.names = FALSE,
                   quote = FALSE, append = TRUE, sep = "\t")
        
        worker_debug_log(paste("Query", query_num, "written to", output_file), "DEBUG")
        
      }, error = function(e) {
        worker_debug_log(paste("ERROR post-processing query", query_num, ":", e$message), "ERROR")
      })
    } else {
      worker_debug_log(paste("Query", query_num, "has no results to write"), "DEBUG")
    }
    
    worker_debug_log(paste("Query", query_num, "completed successfully"), "INFO")
    return(df2)
    
  }, error = function(e) {
    worker_debug_log(paste("FATAL ERROR in query", query_num, ":", e$message), "ERROR")
    return(NULL)
  })
}

# Enhanced parallel processing with proper variable export
if(TRUE){
  no_of_cores <- 30 
  debug_log(paste("Initializing parallel processing with", no_of_cores, "cores"), "INFO")
  
  cl <- makePSOCKcluster(no_of_cores)
  
  total_queries <- length(rownum_query)
  batch_size <- max(1, floor(total_queries / 100))
  
  debug_log(paste("Total queries:", total_queries, "Batch size:", batch_size), "INFO")
  
  # Enhanced cluster setup - define everything workers need
  clusterEvalQ(cl, {
    suppressMessages(library(dplyr))
    options(scipen=999)
  })
  
  debug_log("Exporting variables to cluster...", "INFO")
  
  # Export ALL required variables and functions
  clusterExport(cl, c("rownum_query", "result", "risc", "worker_debug_log"), 
                envir = environment())
  
  # Process batches with detailed debugging
  for(batch_num in 1:100) {
    batch_debug_prefix <- paste("Batch", batch_num)
    debug_log(paste("Starting", batch_debug_prefix), "INFO")
    
    tryCatch({
      start <- (batch_num - 1) * batch_size + 1
      end <- ifelse(batch_num == 100, total_queries, batch_num * batch_size)
      
      debug_log(paste(batch_debug_prefix, "processing queries", start, "to", end), "INFO")
      
      if (start > total_queries) {
        debug_log(paste(batch_debug_prefix, "skipped - start > total_queries"), "INFO")
        next
      }
      
      # Process batch
      results <- parLapply(cl, start:end, ParseEachQuery)
      
      gc()
      debug_log(paste(batch_debug_prefix, "completed successfully"), "INFO")
      
    }, error = function(e) {
      debug_log(paste("ERROR in", batch_debug_prefix, ":", e$message), "ERROR")
      # Try to get more detailed error information
      tryCatch({
        errors <- clusterEvalQ(cl, {
          # Check if there are any recent errors
          if (exists("last_error")) last_error else "No specific error captured"
        })
        debug_log(paste("Worker errors:", paste(unlist(errors), collapse = "; ")), "ERROR")
      }, error = function(e2) {
        debug_log("Could not retrieve worker error details", "ERROR")
      })
    })
  }
  
  stopCluster(cl)
  debug_log("All parallel processing completed", "INFO")
}

debug_log("Script finished", "INFO")