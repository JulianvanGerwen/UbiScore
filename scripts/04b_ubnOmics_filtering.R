# This just stores a few functions for filtering omics data for regulation

#Function: Pivot matrix to long, add metadata, and make psite groups
matrixToLongAndGrouped <- function(matrix){
  #Pivot to long
  matrixLong <- matrix %>% as.data.frame 
  pivCols <- colnames(matrixLong)
  matrixLong <- matrixLong %>%
    rownames_to_column("uniprot_pos") %>%
    pivot_longer(cols = pivCols, names_to = "expCond_id", values_to = "log2FC") %>%
    subset(!is.na(log2FC))
  #Add metadata
  matrixLong <- left_join(matrixLong, ubnOmicsComp_hsPRIDE_lst$meta, by = "expCond_id") %>%
    .[order(.$group), ] %>% 
    mutate(expCond_id = factor(expCond_id, levels = unique(expCond_id))) %>%
    mutate(combinedDescUID = factor(combinedDescUID, levels = unique(combinedDescUID))) %>%
    left_join(., ubnOmicsComp_hsPRIDE_lst$identif, by = "uniprot_pos")
  #Add groups
  matrixLongGrouped <- matrixLong %>%
  dplyr::rename("groupCond" = "group") %>%
  makeCombDat(data = ., args = list(
  "notReg" = list("PSP_reg", F),
  "reg" = list("PSP_reg", T),
  "protDeg" = list("PSP_reg_ProtDeg", T),
  "notProtDeg" = list("PSP_reg_notProtDeg", T)
))
  return(matrixLongGrouped)
}

#Function: classify up/down/either regulated sites by top/bottom fracCutoff quantiles,
#then tally numQuant/numReg/fracReg per condition group into identif
classifyRegulatedSites <- function(ubnOmicsComp_hsPRIDE_lst, fracCutoff = 0.05){
  mat <- ubnOmicsComp_hsPRIDE_lst$matrix
  meta <- ubnOmicsComp_hsPRIDE_lst$meta
  identif <- ubnOmicsComp_hsPRIDE_lst$identif

  #Classify regulated entries using top/bottom fracCutoff
  lower_cutoff <- apply(mat, 2, quantile, probs = fracCutoff, na.rm = TRUE)
  upper_cutoff <- apply(mat, 2, quantile, probs = 1 - fracCutoff, na.rm = TRUE)
  regBoolMLst <- list(
    "down" = sweep(mat, 2, lower_cutoff, FUN = "<"),
    "up" = sweep(mat, 2, upper_cutoff, FUN = ">")
  )
  regBoolMLst$either <- regBoolMLst$down | regBoolMLst$up
  regBoolMLst <- map(regBoolMLst, ~{.x[is.na(.x)] <- F; .x})

  #Conditions to summarise over
  condsLst <- list(
    "all" = meta$expCond_id,
    "notPrtsm" = subset(meta, group != "Proteasome inhibition")$expCond_id,
    "Bort" = subset(meta, group == "Proteasome inhibition" &
                      grepl("Bortezomib", combinedDesc))$expCond_id,
    "BortMG" = subset(meta, group == "Proteasome inhibition" &
                      (grepl("Bortezomib", description) |
                         grepl("MG-132", description)))$expCond_id
  )
  for (grp in unique(meta$group)){
    condsLst[[grp]] <- subset(meta, group == grp)$expCond_id
  }

  #Get # regulated and fraction regulated per condition group and direction
  for (name in names(condsLst)){
    conds <- condsLst[[name]]
    quantNumCol <- paste0("numQuant_", name)
    identif[, quantNumCol] <- rowSums(!is.na(mat[identif$uniprot_pos, intersect(colnames(mat), conds)]))

    for (dir in names(regBoolMLst)){
      suffix <- ifelse(dir == "either", "", paste0("_", dir))
      boolM <- regBoolMLst[[dir]]
      regNumCol <- paste0("numReg_", name, suffix)
      regFracCol <- paste0("fracReg_", name, suffix)

      identif[, regNumCol] <- rowSums(boolM[identif$uniprot_pos, intersect(colnames(boolM), conds)])
      identif[, regFracCol] <- identif[, regNumCol] / identif[, quantNumCol] * 100 #Make percentage
      identif[is.na(identif[, regFracCol]), regFracCol] <- 0
    }
  }

  ubnOmicsComp_hsPRIDE_lst$identif <- identif
  return(ubnOmicsComp_hsPRIDE_lst)
}