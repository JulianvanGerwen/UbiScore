#####Background#####
#Here I store R functions for preparing data for hotspots

library(httr)
library(jsonlite)
library(dplyr)

#####Extracting domains#####
#Function: Extract pfam domains from an Interpro API call
processInterproAPI <- function(data){
  #Identify pfam domains
  #Pfam results seem to just be ordered
  if ("pfam" %in% colnames(data$results$metadata$member_databases)){
    pfamDat <- data$results$metadata$member_databases$pfam
    pfamInds <- which(rowSums(!is.na(pfamDat)) > 0)
    ##Get pfam names
    boolM <- !is.na(pfamDat)
    pfamNames <- sapply(1:nrow(boolM), function(i){paste(colnames(pfamDat)[boolM[i, ]], collapse = ";")})
    pfamNames <- pfamNames[pfamInds]
    #Make output dataframe
    resDf <- data.frame("accession" = data$results$metadata$accession[pfamInds], "pfam" = pfamNames,
                        "source_database" = data$results$metadata$source_database[pfamInds])
    
    #Get positions
    #This accounts for when there are multiple positions listed, by separating the start and end points with ;
    resDf$start <- map(pfamInds, function(i){
      locs <- map(data$results$proteins[[i]]$entry_protein_locations[[1]]$fragments, ~.[1, "start"]) %>% unlist
      return(paste(locs, collapse = ";"))
    }) %>% unlist
    resDf$end <- map(pfamInds, function(i){
      locs <- map(data$results$proteins[[i]]$entry_protein_locations[[1]]$fragments, ~.[1, "end"]) %>% unlist
      return(paste(locs, collapse = ";"))
    }) %>% unlist
    #resDf[, c("start", "end")] <- map(pfamInds, 
    #                                  ~as.numeric(data$results$proteins[[.]]$entry_protein_locations[[1]]$fragments[[1]][1, c("start", "end")])) %>%
    #  purrr::reduce(rbind)
    #NEED TO CHECK IF DISCONTINUOUS
    resDf$dcStatus <- map(pfamInds,
                          ~as.character(data$results$proteins[[.]]$entry_protein_locations[[1]]$fragments[[1]][1, c("dc-status")])) %>% unlist
    return(resDf)
  } else {
    return(NULL)
  }
}


#Function to fetch InterPro entries for a given UniProt protein ID
callInterproAPI <- function(uniprot_id) {
  base_url <- "https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot"
  
  # Construct the API URL
  url <- paste0(base_url, "/", uniprot_id, "?format=json")
  
  # Make the API request
  response <- GET(url)
  
  # Check if the request was successful
  if (status_code(response) != 200) {
    print("Failed to retrieve data. Check UniProt ID and try again.")
    return(NULL)
  } else {
    # Parse JSON response
    data <- fromJSON(httr::content(response, as = "text"), flatten = F)
    return(data)
  }
}


#Function to call and process interpro API
callProcInterproAPI <- function(uniprot_id){
  data <- callInterproAPI(uniprot_id)
  resDf <- processInterproAPI(data = data)
  return(resDf)
}


#####Handling alignments#####
# Function: For a fasta file of alignments, make a dataframe mapping absolute positions to alignment positions for each sequence
# Sequences needed to be named as uniprot;startpos;endpos
get_alignmentMap <- function(fastaF){
  alignment <- read.fasta(fastaF)
  # For each aligned sequence, create a map between absolute position (no gaps) 
  # and alignment position (with gaps)
  alignment_position_maps <- map(names(alignment), function(seq_name) {
    seq_chars <- as.character(alignment[[seq_name]])
    # Get starting positions
    firstPos <- as.numeric(strsplit(seq_name, ";")[[1]][2])
    # Track absolute position (ignoring gaps) and alignment position (including gaps)
    abs_pos <- 0
    mapping_list <- list()
    
    for (aln_pos in 1:length(seq_chars)) {
      if (seq_chars[aln_pos] != "-") {
        abs_pos <- abs_pos + 1
        mapping_list[[abs_pos]] <- aln_pos
      }
    }
    
    # Convert to data frame for easier use
    mapping_df <- data.frame(
      absolute_position = 1:length(mapping_list),
      alignment_position = unlist(mapping_list)
    )
    mapping_df$absolute_position <- mapping_df$absolute_position + firstPos - 1
    
    return(mapping_df)
  })
  names(alignment_position_maps) <- names(alignment)
  return(alignment_position_maps)
}


#####Wrangling results#####
#Function: Give a dataframe a hotspot identifier that uses the hotspot name
get_hotspot_idwName <- function(data, interpro_name = "interpro_name", hotspot_id = "hotspot_id"){
  data <- as.data.frame(data)
  data$hotspot_idwName <- map2(data[, interpro_name], data[, hotspot_id], function(x, y){
    site <- gsub("^.+_", "", y)
    return(paste(x, y, sep = "_"))
  }) %>% unlist 
  return(data)
}









