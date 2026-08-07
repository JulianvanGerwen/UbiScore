# This is code to calculate number of PTMs adjcanet to a ubi-site in PhosphositePlus

library(tidyverse)
library(purrr)
# Load ubiFeatProc
home_directory <- ""
load("data/features/ubiFeatProc__1.RData")

# Load PSP data for all modifications
PSPDir <- paste(home_directory, "data/biol_databases/phosphosite_plus/DATE", sep = "")
Fs <- dir(PSPDir) %>% .[grep("_dataset$", .)]
PSPCombDat <- map(Fs, function(f){
  PSP_raw <- read_delim(file.path(PSPDir, f)) %>%
    as.data.frame
  PSP_sub <- subset(PSP_raw, ORGANISM == "human" & ACC_ID %in% ubiFeatProc$uniprot) %>%
    mutate(mod = gsub("^.+-", "", MOD_RSD),
           site = gsub("-.+$", "", MOD_RSD)) %>%
    mutate(pos = as.numeric(gsub("[A-Z]", "", site))) %>%
    mutate(uniprot_pos = paste(ACC_ID, pos, sep = "_"),
           uniprot_site = paste(ACC_ID, site, sep = "_")) %>%
    dplyr::select("ACC_ID", "mod", "site", "pos", "uniprot_pos", "uniprot_site")
  return(PSP_sub)
}) %>% purrr::reduce(rbind)

# Function to count adjacent PTMs in a window
count_adjacent_PTMs <- function(feat_data, PSP_data, half_window) {
  # Pre-build lookup: for each uniprot ID, store sorted vector of positions
  PSP_by_uniprot <- split(PSP_data$pos, PSP_data$ACC_ID)
  
  # Initialize result
  result <- feat_data %>%
    mutate(pos = as.numeric(pos)) %>%
    dplyr::select("uniprot", "pos", "uniprot_pos") %>%
    mutate(PSP_numAdj = 0)
  
  # Process each row with vectorized range checks
  for (i in 1:nrow(result)) {
    unip <- result$uniprot[i]
    curr_pos <- result$pos[i]
    
    if (unip %in% names(PSP_by_uniprot)) {
      PSP_positions <- PSP_by_uniprot[[unip]]
      # Count positions within window, excluding the position itself
      numAdj <- sum(PSP_positions >= (curr_pos - half_window) & 
                    PSP_positions <= (curr_pos + half_window) & 
                    PSP_positions != curr_pos)
      result$PSP_numAdj[i] <- numAdj
    }
  }
  
  return(result)
}


# Run for multiple window sizes
window_halves <- c(2, 4, 6, 8, 10)
results_by_window <- list()

for (half_w in window_halves) {
  print(paste("Processing window half-width:", half_w))
  start_time <- Sys.time()
  adjPTMRes <- count_adjacent_PTMs(ubiFeatProc, PSPCombDat, half_window = half_w)
  elapsed <- Sys.time() - start_time
  
  col_name <- paste0("PSP_numAdj_", 2*half_w + 1)
  colnames(adjPTMRes)[colnames(adjPTMRes) == "PSP_numAdj"] <- col_name
  results_by_window[[as.character(half_w)]] <- adjPTMRes[, c("uniprot_pos", col_name)]
  
  print(paste("  Time elapsed:", round(elapsed, 2), "seconds"))
  print(paste("  Summary stats for", col_name, ":"))
  print(summary(adjPTMRes[[col_name]]))
}

# Combine
ubiFeatProc_annotated <- ubiFeatProc
for (half_w in window_halves) {
  col_name <- paste0("PSP_numAdj_", 2*half_w + 1)
  ubiFeatProc_annotated <- left_join(ubiFeatProc_annotated, 
                                     results_by_window[[as.character(half_w)]], 
                                     by = "uniprot_pos")
}

windCols <- colnames(ubiFeatProc_annotated) %>% .[grep("^PSP_numAdj_", .)]
PSPwindRes <- ubiFeatProc_annotated[, c("uniprot_site", windCols)]
save(PSPwindRes, file = "output/data/features/PSPwindRes.RData")