###Background
#This script contains functions to assist the use of the R machine learning package "caret"

###Initialise
library(caret)
library(pROC)
library(yardstick)
library(PRROC)
library(ggplot2)
library(tidyverse)
library(purrr)

#####Model metrics and summaries#####

#Function: Extract predictions across CVs for the best version of a model across a hyperparameter search
#instead of supplying a model, you can directly supply prediction data. This is a hack so that ROC plotting code can be used with predictions
extract_CVpreds_bestTune <- function(mod){
  ##Bypass if data is supplied directly
  if ("train" %in% class(mod) == F){
    return(mod)
  ##Bypass if there were no parameter
  } else if (mod$bestTune[1, 1] == "none"){
    return(mod$pred)
  } else {
    ##Get out predictions for all parameters and CVs
    predAll <- mod$pred
    predBest <- predAll
    ##Subset for the best model
    for (col in colnames(mod$bestTune)){
      val <- mod$bestTune[1, col]
      predBest <- predBest[which(predBest[, col] == val), ]
    }
    return(predBest)
  }
}


#####ROC curves and related#####

####ROC curves using yardstick
#I use the packages yardstick (ROC) or PRROC (PR) to create dataframes that can be used to make plots
#For multiple folds I extrapolate these dataframes so that they share the same x values. Then I combine

###Function: Makes data for plotting a ROC curve or other model curve given a model with cross validation. 
#Will interpolate individual curves to get mean and sd
#mod: The model
#type: ROC or PR
#output: if MeanSd, return wide data with mean and sd at each x value
#if long, return long data with folds indicated in a column Resample
prep_ROC_meansd <- function(mod,
                            type = "ROC",
                            output = "MeanSd"){
  ##Get curve data for each fold
  #This depends on the type of curve
  if (type == "ROC"){
    curve_fold_lst <- extract_CVpreds_bestTune(mod) %>% split(.$Resample) %>%
      map2(., names(.), function(dat, fold){
        ROC_dat <- roc_curve(data = dat, truth = "obs", "X1")
        ROC_dat$OneMinSpec <- 1 - ROC_dat$specificity
        #Rename x and y for generality
        ROC_dat <- dplyr::rename(ROC_dat,
                                 "x" = "OneMinSpec",
                                 "y" = "sensitivity")
        return(ROC_dat)
      })
  } else if (type == "PR"){ #I use the package PRROC because it uses nonlinear interpolation
    curve_fold_lst <- extract_CVpreds_bestTune(mod) %>% split(.$Resample) %>%
      map2(., names(.), function(dat, fold){
        #Prepare data
        dat$obs_num <- (dat$obs == "X1")*1
        #Make pr object
        pr <- pr.curve(scores.class0 = dat$X1,
                       weights.class0 = dat$obs_num,
                       curve = TRUE)
        #Extract values
        pr_df <- pr$curve %>% as.data.frame
        colnames(pr_df) <- c("x", "y", "colour")
        pr_df <- pr_df[, c("x", "y")]
        return(pr_df)
  })
  }
  
  #Make long data for output
  curve_fold_long <- map2(curve_fold_lst, names(curve_fold_lst), function(dat, fold){
    dat$Resample <- fold
    return(dat)
  }) %>% purrr::reduce(rbind) %>%
    mutate(Resample = factor(Resample, levels = unique(Resample)))
  
  ##Interpolate across all x values
  x_all <- map(curve_fold_lst, ~.$x) %>% unlist %>% unique %>% sort
  ##Interpolate
  curve_fold_lst_int <- map(curve_fold_lst, function(dat){
    int <- approx(dat$x, dat$y, 
                  xout = x_all, ties = "max") %>% #When there is a tie I use the max y-value
      purrr::reduce(cbind) %>% as.data.frame
    colnames(int) <- c("x", "y")
    return(int)
  })
  ##Combine
  curve_fold_int <- map2(curve_fold_lst_int, names(curve_fold_lst_int), function(dat, fold){
    dat$Resample <- fold
    return(dat)
  }) %>% purrr::reduce(rbind) %>%
    mutate(Resample = factor(Resample, levels = unique(Resample)))
  ##Get mean and SD
  folds <- levels(curve_fold_int$Resample)
  curve_fold_int_wide <- pivot_wider(curve_fold_int, 
                                     id_cols = "x",
                                     names_from = "Resample",
                                     values_from = "y")
  curve_fold_int_wide$mean <- apply(curve_fold_int_wide[, folds], 1, "mean", na.rm = TRUE)
  curve_fold_int_wide$sd <- apply(curve_fold_int_wide[, folds], 1, "sd", na.rm = TRUE)
  ##Add back the origin if it is missing
  if (type == "ROC"){
    if (0 %in% curve_fold_int_wide$mean == FALSE){
      curve_fold_int_wide <- rbind(rep(0, nrow(curve_fold_int_wide)), curve_fold_int_wide)
    } 
  }
  
  ##Output
  if (output == "MeanSd"){
    return(curve_fold_int_wide)
  } else if (output == "long"){
    return(curve_fold_long)
  }
}

##Prep data for multiple models. See above for arguments
#mods must be a named list of models
prep_ROC_meansd_fromLst <- function(mods,
                                    type = "ROC",
                                    output = "MeanSd"){
  ##Account for if only one model is supplied
  if (class(mods)[1] != "list"){
    mods <- list("model" = mods)
  }
  ##Prep data
  dataLst <- mods %>% map2(., names(.), function(mod, name){
    dat <- prep_ROC_meansd(mod = mod, type = type, output = output)
    dat$model <- name
    return(dat)
  }) 
  ##Pad with extra columns in case they are missing. This happens if a resample has been thrown out
  commonCols <- map(dataLst, colnames) %>% unlist %>% unique
  dataLst <- map(dataLst, function(dat){
    missingCols <- setdiff(commonCols, colnames(dat))
    if (length(missingCols) > 0){
      for (col in missingCols){
        dat[, col] <- NA
      }
    }
    dat <- dat[, commonCols]
    return(dat)
  })
  ##Combine
  data <- dataLst %>% purrr::reduce(rbind) %>%
    mutate(model = factor(model, levels = unique(model))) 
  return(data)
}


###Function: Plot ROC curve given data output by prep_ROC_meansd_fromLst
#type: ROC or PR. Determines axis labels
##AUCdat: data with columns "AUC" and the same columns as "fillAes" and "colourAes", as factors. AUC values are added to top left corners
plotROC_fromData <- function(data,
                             colourAes = "model", fillAes = "model",
                             ribbonAlpha = 0.15,
                             MeanSd = T,
                             type = "ROC", coordFixed = T,
                             AUCdat = NULL,
                             linewidth = 0.235){
  ##Set up for mean +-sd
  if (MeanSd){
    data$y <- data$mean
    ribbon_plot <- geom_ribbon(alpha = ribbonAlpha,
                               colour = NA)
    path_plot <- geom_path(linewidth = linewidth)
  } else {
    if ("y" %in% colnames(data) == F){ #Account for when data is prepped for meansd, but we don't want SD
      data$y <- data$mean
    }
    data$sd <- 0
    data$groupAes <- paste(data$Resample, data$model, sep = "_")
    ribbon_plot <- NULL
    path_plot <- geom_path(linewidth = linewidth,
                           aes(group = groupAes))
  }
  ##Rename aesthetics
  renames <- c("fillAes" = fillAes, "colourAes" = colourAes)
  for (i in 1:length(renames)){
    data[, names(renames[i])] <- data[, renames[i]]
  }
  ##Set up abline
  #Set up reduced data for abline
  #This contains every observed combination of levels for columns that are factors. 
  #Hence, the abline is plotted once for each combination of levels
  factorColBool <- map(colnames(data), ~is.factor(as.data.frame(data)[, .])) %>% unlist
  factorCols <- colnames(data)[factorColBool]
  data$factorsComb <- apply(data[, factorCols], 1, function(x){paste(x, collapse = "_")})
  datSum <- data %>% split(.$factorsComb) %>% map(~.[1, ]) %>% purrr::reduce(rbind)
  if (type == "ROC"){
    abline_plot <- geom_segment(data = datSum,
                                x = 0, y = 0, xend = 1, yend = 1,
                                colour = "black", alpha = 0.7,
                                linewidth = 0.235, linetype = "dashed")
  } else {
    abline_plot <- NULL
  }
  ##Set up labels
  if (type == "ROC"){
    labels <- labs(x = "1 - Specificity",
                   y = "Sensitivity")
  } else if (type == "PR"){
    labels <- labs(x = "Recall", y = "Precision")
  }
  ##Set up scales
  if (coordFixed == T){
    coords <- coord_fixed()
  } else {
    coords <- NULL
  }
  ##Set up AUC text
  if (!is.null(AUCdat)){
    ##Set up AUC plot
    AUCdatPlot <- AUCdat %>% mutate(x = 0.05, y = 0.95, sd = 0,
                                    AUC = signif(AUC, 3))
    #If fillAes is not present, add it - otherwise colours break
    if (fillAes %in% colnames(AUCdatPlot) == F){
      AUCdatPlot[, fillAes] <- as.character(as.data.frame(ROCdat)[1, fillAes])
      AUCdatPlot[, fillAes] <- factor(AUCdatPlot[, fillAes], levels = levels(as.data.frame(ROCdat)[, fillAes]))
    }
    #Make text plot
    textPlot <- geom_text(data = AUCdatPlot, aes(label = AUC, x = x, y = y, fill = model),
                          size = 1.875, hjust = 0, colour = "black")
  } else {
    textPlot <- NULL
  }
  ##Plot
  output_plot <- ggplot(data, aes(x = x, y = y,
                                  ymax = y + sd, ymin = y - sd,
                                  colour = colourAes, fill = fillAes)) +
    ribbon_plot +
    path_plot + 
    abline_plot +
    textPlot +
    comfy_theme() +
    coords +
    labels
  return(output_plot)
}

###Function: Plot ROC curves from a list of models
plotROC_fromMods <- function(mods, type = "ROC",
                             colourAes = "model", fillAes = "model",
                             ribbonAlpha = 0.15,
                             MeanSd = T,
                             ...){
  ##Prep data
  if (MeanSd){
    output <- "MeanSd"
  } else {
    output <- "long"
  }
  data <- prep_ROC_meansd_fromLst(mods = mods,
                                  type = type, output = output)
  ##Plot
  output_plot <- plotROC_fromData(data = data,
                                  colourAes = colourAes, fillAes = fillAes,
                                  ribbonAlpha = ribbonAlpha,
                                  MeanSd = MeanSd,
                                  type = type,
                                  ...)
  return(output_plot)
}


###Function: Plot ROC curve from data with a binary outcome and a numerical predictor column
#obsCol: Column referring to the binary outcome. Factor with levels X0 and X1, where X1 is positive case
#predCol: A numerical col predicting the binary outcome. Higher values try to predict X1
#groupAes: A column to use for grouping multiple ROC curves that will be overlaid
plotROC_fromPredictor <- function(data, obsCol, predCol, colour = "black", linewidth = 0.235,
                                  groupAes = NULL){
  ##Prep dat
  ###Use group aesthetic for multiple curves
  data[, c("obs", "pred")] <- data[, c(obsCol, predCol)]
  if (is.null(groupAes)){
    data$group <- 1
  } else {
    data$group <- data[, groupAes]
  }
  plotDat <- data %>% split(.$group) %>% map(~roc_curve(data = ., truth = "obs", "pred"))
  plotDat <- map2(plotDat, names(plotDat), function(dat, group){
    dat$group <- group
    return(dat)
  }) %>% purrr::reduce(rbind) %>%
    mutate(group = factor(group, levels = unique(group)))
  
  ##Set up abline
  #Set up reduced data for abline
  #This contains every observed combination of levels for columns that are factors. 
  #Hence, the abline is plotted once for each combination of levels
  factorColBool <- map(colnames(plotDat), ~is.factor(as.data.frame(plotDat)[, .])) %>% unlist
  factorCols <- colnames(plotDat)[factorColBool]
  plotDat$factorsComb <- apply(plotDat[, factorCols], 1, function(x){paste(x, collapse = "_")})
  datSum <- plotDat %>% split(.$factorsComb) %>% map(~.[1, ]) %>% purrr::reduce(rbind)
  abline_plot <- geom_segment(data = datSum,
                                  x = 0, y = 0, xend = 1, yend = 1,
                                  colour = "black", alpha = 0.7,
                                  linewidth = 0.235, linetype = "dashed")
  ##Use group aesthetic to overlay curves
  if (!is.null(groupAes)){
    output_plot <- ggplot(plotDat, aes(x = 1 - specificity, y = sensitivity, group = group, colour = group)) +
      geom_path(linewidth = linewidth) +
      abline_plot +
      comfy_theme() +
      coord_fixed()
  } else {
    output_plot <- ggplot(plotDat, aes(x = 1 - specificity, y = sensitivity)) +
      geom_path(linewidth = linewidth, colour = colour) +
      abline_plot +
      comfy_theme() +
      coord_fixed()
  }
  return(output_plot)
}

###Function: Calculate ROC AUC from data with a binary outcome and a numerical predictor column
#obsCol: Column referring to the binary outcome. Factor with levels X0 and X1, where X1 is positive case
#predCol: A numerical col predicting the binary outcome. Higher values try to predict X1
ROCAUC_fromPredictor <- function(data, obsCol, predCol){
  ##Prep dat
  data[, c("obs", "pred")] <- data[, c(obsCol, predCol)]
  ##Make obsCol a factor if it is logical
  if (typeof(data$obs) == "logical"){
    data$obs <- factor(make.names(as.numeric(data$obs)), levels = c("X1", "X0"))
  }
  ##Get AUC
  AUCDat <- roc_auc(data, truth = "obs", "pred")
  return(as.numeric(AUCDat[1, ".estimate"]))
}

# Calcualte ROC AUC from vectors. These must be formatted like columns that would go into ROCAUC_fromPredictor
ROCAUC_fromVecs <- function(obsVec, predVec){
	data <- data.frame("obs" = obsVec, "pred" = predVec)
	return(ROCAUC_fromPredictor(data = data, obsCol = "obs", predCol = "pred"))
}