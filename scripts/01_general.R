#####Background#####
#Here I store general scripts for manipulating data etc.


#####PhosphositePlus#####
#Split PSP regulatory site data by ON_function, so each row is a unique site-ON_FUNCTION_dir pair
split_ON_FUNCTION <- function(data){
  PSP_reg_ub_FUNC <- data %>% subset(!is.na(ON_FUNCTION)) %>%
    separate_rows(ON_FUNCTION, sep = "; ")
  PSP_reg_ub_FUNC$ON_FUNCTION_strp <- strsplit(PSP_reg_ub_FUNC$ON_FUNCTION, ", ") %>% map(~.[1]) %>% unlist
  PSP_reg_ub_FUNC$ON_FUNCTION_dir <- strsplit(PSP_reg_ub_FUNC$ON_FUNCTION, ", ") %>% map(function(x){
    if (length(x) == 1){return(NA)} else {return(x[2])}
  }) %>% unlist
  return(PSP_reg_ub_FUNC)
}



#####Data summarising#####
#Function: Summarise the number and percentage of true values in boolean columns over a categorical column
catBoolSummariser <- function(data, catCol, boolCols){
  summary_df <- data %>%
    group_by(!!sym(catCol)) %>%
    summarise(
      num_rows = n(),  # Total rows in each group
      across(all_of(boolCols), 
             list(
               num = ~ sum(.),  # Count TRUE values
               pct = ~ sum(.) / n() * 100  # Percentage of TRUE
             ), 
             .names = "{.col}_{.fn}") # Add suffix to column names
    )
  ##Calculate % of the TRUE values in each category value
  for (col in boolCols){
    numCol <- paste0(col, "_num")
    summary_df[, paste0(col, "_groupPct")] <- summary_df[, numCol] / sum(summary_df[, numCol]) * 100
  }
  ##Calculate enrichment of TRUE values in each category value
  for (col in boolCols){
    numCol <- paste0(col, "_num")
    pctCol <- paste0(col, "_pct")
    summary_df[, paste0(col, "_enr")] <- summary_df[, pctCol] / (sum(summary_df[, numCol])/sum(summary_df[, c("num_rows")])) / 100
  }
  return(summary_df)
}


#Function: Paste the same data subsetted in different ways. Good for violin plots
makeCombDat <- function(data, args = list(
  "notReg" = list("PSP_reg", F),
  "protDeg" = list("PSP_reg_ProtDeg", T),
  "otherFunc" = list("PSP_reg_otherFunc", T)
)){
  datLst <- map(args, function(x){
    filter(data, !!sym(x[[1]]) == x[[2]])
  })
  combDat <- map_dfr(datLst, function(x){x}, .id = "group") %>%
    mutate(group = factor(group, levels = names(datLst)))
  return(combDat)
}


#Function: Summarise Group a numerical column into bins
binner <- function(data, bin_breaks = c(0, 5, 10, 20, Inf), valCol,
                   breakType = "integer"){
  # Generate bin labels dynamically
  if (breakType == "integer"){
    bin_labels <- paste0(head(bin_breaks, -1), "-", tail(bin_breaks, -1) - 1)
  } else {
    bin_labels <- paste0(head(bin_breaks, -1), "-", tail(bin_breaks, -1))
  }
  ##If Inf in there
  if (Inf %in% bin_breaks){bin_labels[length(bin_labels)] <- paste0(bin_breaks[length(bin_breaks) - 1], "+")}

  
  # Apply binning and count within groups
  df_binned <- data %>%
    mutate(Bin = cut(!!sym(valCol), breaks = bin_breaks, labels = bin_labels, right = F, include.lowest = TRUE)) 
  return(df_binned)
}


#Function: Group a numerical column into bins and then count number of entries in each bin, in each value of a group column
#The bin breaks denote the starting point of each break. The bins are closest on the left and open on the right. 
#So in the given example, the first two bins are [0, 6), [6, 11)
binAndSummarise <- function(data, bin_breaks = c(0, 6, 11, 21, Inf),
                            groupCol, valCol,
                            breakType = "integer"){
  df_binned <- binner(data = data, 
                      bin_breaks = bin_breaks, 
                      valCol = valCol,
                      breakType = breakType) %>%
    dplyr::count(!!sym(groupCol), Bin) %>%
    group_by(!!sym(groupCol)) %>%
    mutate(Percentage = 100 * n / sum(n)) %>%
    ungroup()
  
  # Print result
  return(df_binned)
}

# Functional: Cut a numerical vector into bins, where the breaks are quantiles specified in quantBreaks
# E.g. quantBreaks = c(0, 0.5, 0.1) will split into two bins separated at the median
cutWithQuantBreaks <- function(x, quantBreaks){
  xCut <- cut(x, quantile(x, quantBreaks, na.rm = T))
  # Account for values at the bottom that are cut off
  xCut[which(x == min(x, na.rm = T))] <- levels(xCut)[1]
  return(xCut)
}



##### Maths #####
geoMean <- function(x){exp(mean(log(x)))}

# Convert vector to quantiles
to_quantiles <- function(x){
  if (sum(!is.na(x)) == 0){
    return(x)
  } else {
    return(ecdf(x)(x))
  }
} 
