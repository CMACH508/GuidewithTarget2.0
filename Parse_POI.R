args=commandArgs(T)
file1 = args[1] # POI_Target_Pair:ago3_POI_Target_Pair.txt
#data_add = args[2]
#reference = args[2]
risc = as.character(args[2] )
#reference = ifelse(reference=="1","transcriptome","genome")
options(scipen=999)

#print(file1)
#print(reference)
#print(iteration)
#print(risc)


suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library(parallel))




#file1 = "ago2_guide1_POI_BP.txt"
#data_add = "/data3/zenglin/temp/LinZeng/ZhongJing_degradome/wkdir/ago2_guide1_Cleavage1_NoCleavage1/blast_POI/"
#reference = "2"
#risc = "ago2_guide1"

# Read the file with error handling
result <- tryCatch({
  readLines(paste0(file1))
}, error = function(e) {
  stop(paste("Error reading file:", e$message))
})
rownum_query <- which(startsWith(result,"Query="))


#system(paste0("rm ",data_add,risc,"_Parsed_NPOI_",iteration,".txt"))
#system(paste0("rm ", data_add,risc,"_Parsed_Filter_NPOI_",iteration,".txt"))

# Enhanced ParseEachQuery function with better error handling
ParseEachQuery <- function(query_num){
  tryCatch({
    # Check if query_num is valid
    if (query_num < 1 || query_num > length(rownum_query)) {
      return(data.frame())  # Return empty dataframe for invalid query_num
    }
    
    removesym <- function(i){
      pro_i <- substr(i,3,nchar(i))
      return(pro_i)
    }
    
    # Your existing function code here, but with additional safety checks
    df2 <- data.frame(query_id=character(), sub_id=character(), q_start=character(),
                      q_end=character(), s_start=character(),  s_end=character(),
                      q_match_pos=character(), q_mismatch_pos=character(),
                      identity=character(), alignment_length = numeric(),  match_num = numeric(),
                      mismatch_num = numeric(), stringsAsFactors=FALSE)
    
    query_row = rownum_query[query_num]
    
    # Safer extraction of query_extrat
    if(length(rownum_query) == 1){
      query_extrat <- result[query_row:length(result)]
    } else if(query_num == length(rownum_query)){
      query_extrat <- result[query_row:length(result)]
    } else {
      next_query_start <- rownum_query[query_num + 1]
      # Ensure we don't go beyond result length
      if (next_query_start - 1 <= length(result)) {
        query_extrat <- result[query_row:(next_query_start - 1)]
      } else {
        query_extrat <- result[query_row:length(result)]
      }
    }
    
    # Check if query_extrat is valid
    if (length(query_extrat) == 0) {
      return(data.frame())
    }
    
    sub_id_list = query_extrat[which(startsWith(query_extrat,">"))]
    query_id = substr(query_extrat[1],8,nchar(query_extrat[1]))
    Nohist = query_extrat[6]
    
    # Check if Nohist exists
    if (length(Nohist) == 0) {
      return(data.frame())
    }
    
    if(Nohist=="***** No hits found *****"){
      print(paste("No Hits for query", query_num))
    } else {
      # Check if sub_id_list is not empty
      if (length(sub_id_list) == 0) {
        return(data.frame())
      }
      
      for(subject in 1:length(sub_id_list)){
        sub_id = sub_id_list[subject]
        sub_rownum = which(startsWith(query_extrat, sub_id))
        
        # Check if sub_rownum exists and is valid
        if (length(sub_rownum) == 0 || sub_rownum > length(query_extrat)) {
          next
        }
        
        sub_rownum1 = min(sub_rownum + 11, length(query_extrat))  # Ensure we don't exceed bounds
        content_rel <- query_extrat[sub_rownum:sub_rownum1]
        
        # Check if we have enough content
        if (length(content_rel) < 5) {
          next
        }
        
        identity <- tryCatch({
          (unlist(strsplit(content_rel[5],"), Gaps = "))[1] %>% strsplit("\\(") %>% unlist)[2]
        }, error = function(e) NA)
        
        alignment_summary = tryCatch({
          (unlist(strsplit(content_rel[5]," Identities = "))[2] %>% strsplit(" \\(") %>% unlist)[1]
        }, error = function(e) NA)
        
        alignment_length = tryCatch(as.numeric(unlist(strsplit(alignment_summary,"/"))[2]), error = function(e) NA)
        match_num = tryCatch(as.numeric(unlist(strsplit(alignment_summary,"/"))[1]), error = function(e) NA)
        mismatch_num = tryCatch(alignment_length - match_num, error = function(e) NA)
        
        temp1 = unlist(strsplit(content_rel[8]," "))
        temp1 = temp1[which(temp1!="")]
        q_start <- if(length(temp1) >= 4) temp1[2] else NA
        q_end <- if(length(temp1) >= 4) temp1[4] else NA
        
        temp2 = unlist(strsplit(content_rel[10]," "))
        temp2 = temp2[which(temp2!="")]
        s_end <- if(length(temp2) >= 4) temp2[2] else NA
        s_start <- if(length(temp2) >= 4) temp2[4] else NA
        
        matchpos <- content_rel[9]
        q_match_pos = ""
        q_mismatch_pos = ""
        
        if (!is.na(matchpos) && nchar(matchpos) > 0) {
          count_padding = 0
          while(count_padding < nchar(matchpos) && substr(matchpos, count_padding+1, count_padding+1) != "|") {
            count_padding = count_padding + 1
          }
          for(i in 1:nchar(matchpos)){
            symbol <- substr(matchpos,i,i)
            if(symbol=="|"){
              q_match_pos = c(q_match_pos, i-count_padding+as.numeric(q_start)-1)
            }
          }
          q_match_pos = q_match_pos[2:length(q_match_pos)]
          if (!is.na(q_start) && !is.na(q_end)) {
            q_mismatch_pos = setdiff(c(as.numeric(q_start):(as.numeric(q_end)-1)), q_match_pos)
          }
        }
        
        q_match_pos <- paste0(q_match_pos, collapse=",")
        q_mismatch_pos <- paste0(q_mismatch_pos, collapse=",")
        
        pairresult <- data.frame(query_id=query_id, sub_id=sub_id, q_start=q_start, q_end=q_end, 
                                 s_start=s_start, s_end=s_end, q_match_pos=q_match_pos, 
                                 q_mismatch_pos=q_mismatch_pos, alignment_length=alignment_length,
                                 match_num=match_num, mismatch_num=mismatch_num, identity=identity)
        df2 <- rbind(df2, pairresult)
      }
    }
    
    # Your existing processing code for df2...
    if (nrow(df2) > 0) {
      df2$dis1 <- 10 - as.numeric(df2$q_start)
      pattern_before_plus = ".*_(\\d+)_[+/-].*"
      df2$pos1_based <- as.numeric(gsub(pattern_before_plus, "\\1", df2$sub_id))
      pattern_before_parentheses = ".*-(\\d+)\\([+-]\\).*"
      df2$start_abs <- as.numeric(gsub(pattern_before_parentheses,"\\1", df2$sub_id))
      relpos <- abs(df2$pos1_based - df2$start_abs)
      dis <- relpos - as.numeric(df2$s_end)
      df2$rel_start <- as.numeric(df2$s_start) - relpos
      df2$rel_end <- as.numeric(df2$s_end) - relpos
      
      write.table(file=paste0(risc, "_Parsed_POI_",query_num,".txt"), df2, row.names=FALSE, 
                  col.names=FALSE, quote=FALSE, append=TRUE, sep="\t")
    }
    
    print(paste0("Query ", query_num, " finished!"))
    return(df2)
    
  }, error = function(e) {
    print(paste("Error in query", query_num, ":", e$message))
    return(data.frame())  # Return empty dataframe on error
  })
}

# Modified parallel execution with better error handling
if(TRUE){
  no_of_cores <- 16 
  cl <- makePSOCKcluster(no_of_cores)
  
  # Calculate batch size
  total_queries <- length(rownum_query)
  batch_size <- max(1, floor(total_queries / 100))
  
  print(paste0("100 batches, ", total_queries, " rows to parse!"))
  
  # Set up cluster
  clusterEvalQ(cl, {
    suppressMessages(library(dplyr))
  })
  
  clusterExport(cl, c("rownum_query", "result", "risc"))
  
  # Process in smaller batches with better error handling
  for(batch_num in 1:100) {
    print(paste("Processing batch", batch_num))
    start = (batch_num - 1) * batch_size + 1
    end = ifelse(batch_num == 100, total_queries, batch_num * batch_size)
    
    # Skip if start > end
    if (start > total_queries) break
    
    # Use tryCatch for the entire batch
    tryCatch({
      results <- parLapply(cl, start:end, ParseEachQuery)
      gc()
    }, error = function(e) {
      print(paste("Error in batch", batch_num, ":", e$message))
    })
  }
  
  stopCluster(cl)
}





