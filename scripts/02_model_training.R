#####Background#####
#Here I store objects and functions for training models


#####Preparing features and response#####
###Combine features and resp
combine_FeatResp <- function(feats, resp){
  ##Combine features if they are a list
  if (class(feats) == "list"){
    feats <- purrr::reduce(feats, cbind)
  }
  ##Add response
  comb <- cbind(resp, feats)
  colnames(comb)[1] <- "response"
  return(comb)
}

#####Subsampling negative class#####
##### Subsampling negative class#####
###Function: Subsample negative classes to be a multiple of positives
#resp: A response vector output by respPrep_proc
#subsampNum: Number of subsamples to produce
#subsampMult: How many times bigger than positive classes should negatives be?
subsample_negatives <- function(resp,
                                subsampSeed = 66,
                                subsampNum = 3,
                                subsampMult = 2,
                                reorder = TRUE){
  ##Split positive and negative
  pos_ind <- which(resp == "X1")
  neg_ind <- which(resp == "X0")
  ##Sample from negative
  set.seed(subsampSeed)
  neg_samps <- map(1:subsampNum,
                   ~sample(neg_ind, size = subsampMult*length(pos_ind)))
  ##Reconstitute
  samps <- map(neg_samps, function(neg_samp){
    rec <- c(pos_ind, neg_samp)
    if (reorder){ #if desired, scramble indices
      rec <- sample(rec, size = length(rec))
    }
  }) 
  names(samps) <- paste("Subsamp", 1:subsampNum, sep = "")
  return(samps)
}


###Function: Train a model using subsample_negatives
#All arguments from subsample_negatives are inherited
#modFunc must be a "modTrain" function.
#... is used to pass arguments to modFunc
subsample_train <- function(feats, resp, 
                            subsampSeed = 66,
                            subsampNum = 3,
                            subsampMult = 2,
                            reorder = TRUE,
                            modFunc,
                            ...){
  ##Peform subsampling
  subsamps <- subsample_negatives(resp = resp,
                                  subsampSeed = subsampSeed,
                                  subsampNum = subsampNum,
                                  subsampMult = subsampMult,
                                  reorder = reorder)
  
  ##Train models
  mods <- map(subsamps, function(subsamp){
    ##Apply subsample to features and response
    respSubsamp <- resp[subsamp]
    if (class(feats) == "list"){
      featsSubsamp <- map(feats, ~dplyr::slice(., subsamp))
    } else {
      featsSubsamp <- dplyr::slice(feats, subsamp)
    }
    ##Train model
    mod <- modFunc(feats = featsSubsamp,
                   resp = respSubsamp,
                   ...)
    return(mod)
  })
  
  ##Combine results into unified model
  mods_comb <- mods[[1]] #Copy the first model as a template
  #Keep all models
  mods_comb$Model_subsamps <- mods
  #Extract CV measures
  CV_names <- c("pred", "trainingData", "resample", "resampledCM")
  for (name in CV_names){
    datLst <- map(mods, ~.[[name]])
    datLst <- map2(datLst, names(datLst), function(dat, name){
      dat$Resample <- paste(name, dat$Resample, sep = "")
      return(dat)
    })
    datComb <- purrr::reduce(datLst, rbind)
    mods_comb[[name]] <- datComb
  }
  #trainingData - correct rowIndex
  for (i in 1:subsampNum){
    subsamp_tag <- paste("Subsamp", i, "Fold", sep = "")
    mods_comb$pred$rowIndex[grep(subsamp_tag, mods_comb$pred$Resample)] <- mods_comb$pred$rowIndex[grep(subsamp_tag, mods_comb$pred$Resample)] +
      (i - 1)*nrow(mods_comb$Model_subsamps[[1]]$trainingData)
  }
  #Average results
  ##Here I just average across the same hyperparamter values, rather than the best model from each
  dataList <- map(mods, ~.$results)
  # Combine results from each subsample in a column-wise manner;
  # average numeric columns and paste non-numeric columns using ';' if values differ.
  n <- length(dataList)
  combinedData <- dataList[[1]]
  for(colName in names(combinedData)) {
    if(is.numeric(combinedData[[colName]])) {
      # Average numeric columns
      combinedData[[colName]] <- Reduce("+", lapply(dataList, function(df) df[[colName]])) / n
    } else {
      # Combine non-numeric columns row-by-row:
      # If all values are the same, retain it; otherwise, paste unique values separated by ';'
      combinedData[[colName]] <- sapply(seq_len(nrow(combinedData)), function(i) {
        vals <- sapply(dataList, function(df) as.character(df[[colName]][i]))
        if(length(unique(vals)) == 1) {
          unique(vals)[1]
        } else {
          paste(unique(vals), collapse = ";")
        }
      })
    }
  }
  mods_comb$results <- combinedData
  return(mods_comb)
}


#####Training#####
#Class weights
caseWeights_balanced <- function(resp){
  tab <- table(resp)
  weights <- sapply(resp, function(x){
    return(0.5/tab[which(names(tab) == x)])
  })
  return(weights)
}

#Repeated cross validation
fitControl_repCV <- trainControl(
  #5-fold cross validation
  method = "repeatedcv", repeats = 3, number = 5,
  #Necessary for ROC
  classProbs = T,
  #Misc parameters
  savePredictions = T,
  #Set this to use ROC
  summaryFunction = accROCSummary
)

#This is simply five-fold cross validation
fitControl_CV <- trainControl(
  #5-fold cross validation
  method = "cv", number = 5,
  #Necessary for ROC
  classProbs = T,
  #Misc parameters
  savePredictions = T,
  #Set this to use ROC
  summaryFunction = accROCSummary
)

#Function: Make splits with class balance of "response" preserved
makeFolds <- function(seed = 123, k = 5, data){
  set.seed(seed)  # for reproducibility
  # Create 5-fold stratified train-test splits (each with ~80% train and ~20% test)
  folds <- caret::createFolds(data$response, k = k, returnTrain = TRUE)
  trainTestSplits <- lapply(folds, function(trainIndex) {
    list(
      trainDat = data[trainIndex, ],
      testDat  = data[-trainIndex, ]
    )
  })
  names(trainTestSplits) <- paste0("Fold", 1:length(trainTestSplits))
  return(trainTestSplits)
}

#Function: Make splits with class balance of "response" preserved, stratified across different proteins
makeFoldsProtstrat <- function(seed = 123, k = 5, data){
  ## Add uniprot
  data$uniprot <- gsub("_.+$", "", rownames(data))

  #Collapse at the protein level so I can stratify proteins
  respMax <- function(x){
    if ("X1" %in% x){return(T)} else {return(F)}
  }
  trainDat_protAgg <- data[, c("uniprot", "response")] %>% 
                    group_by(uniprot) %>% summarise(response = respMax(response))

  #Make splits on protein level 
  set.seed(seed)  # for reproducibility
  # Create 5-fold stratified train-test splits (each with ~80% train and ~20% test)
  folds_protein <- caret::createFolds(trainDat_protAgg$response, k = k, returnTrain = TRUE)
  trainTestSplits_protein <- lapply(folds_protein, function(trainIndex) {
    list(
      trainDat = trainDat_protAgg[trainIndex, ],
      testDat  = trainDat_protAgg[-trainIndex, ]
    )
  })
  names(trainTestSplits_protein) <- paste0("Fold", 1:length(trainTestSplits_protein))
  # get site-level folds
  trainTestSplits <- lapply(trainTestSplits_protein, function(x){
    list(
      trainDat = subset(data, uniprot %in% x$trainDat$uniprot) %>% dplyr::select(-c("uniprot")),
      testDat = subset(data, uniprot %in% x$testDat$uniprot) %>% dplyr::select(-c("uniprot"))
    )
  })
  return(trainTestSplits)
}


#####Models#####
modTrain_LogReg_weighted <- function(feats, resp, 
                            fitControl, seed,
                            weightFunc = caseWeights_balanced){
  ##Generate weights
  weights <- weightFunc(resp)
  ##Combine features and resp
  FeatResp_data <- combine_FeatResp(feats = feats, resp = resp)
  ##Train model
  set.seed(seed)
  mod <- train(response ~ ., data = FeatResp_data,
               method = "glm", family = "binomial",
               metric = "Accuracy", trControl = fitControl,
               weights = weights)
  return(mod)
}

modTrain_GBM_weighted <- function(feats, resp, 
                            fitControl, seed,
                            weightFunc = caseWeights_balanced){
  ##Generate weights
  weights <- weightFunc(resp)
  ##Combine features and resp
  FeatResp_data <- combine_FeatResp(feats = feats, resp = resp)
  ##Train model
  set.seed(seed)
  mod <- train(response ~ ., data = FeatResp_data,
                method = "gbm",
                metric = "ROC", trControl = fitControl,
                weights = weights, verbose = F)
  return(mod)
}

modTrain_RF_weighted <- function(feats, resp, 
                            fitControl, seed,
                            weightFunc = caseWeights_balanced){
  ##Generate weights
  weights <- weightFunc(resp)
  ##Combine features and resp
  FeatResp_data <- combine_FeatResp(feats = feats, resp = resp)
  ##Train model
  set.seed(seed)
  mod <- train(response ~ ., data = FeatResp_data,
                method = "ranger",
                metric = "ROC", trControl = fitControl,
                weights = weights)
  return(mod)
}

modTrain_LogReg_unweighted <- function(feats, resp, 
                            fitControl, seed){
  ##Combine features and resp
  FeatResp_data <- combine_FeatResp(feats = feats, resp = resp)
  ##Train model
  set.seed(seed)
  mod <- train(response ~ ., data = FeatResp_data,
               method = "glm", family = "binomial",
               metric = "Accuracy", trControl = fitControl)
  return(mod)
}

modTrain_GBM_unweighted <- function(feats, resp, 
                            fitControl, seed){
  ##Combine features and resp
  FeatResp_data <- combine_FeatResp(feats = feats, resp = resp)
  ##Train model
  set.seed(seed)
  mod <- train(response ~ ., data = FeatResp_data,
                method = "gbm",
                metric = "ROC", trControl = fitControl,
                verbose = F)
  return(mod)
}

modTrain_RF_unweighted <- function(feats, resp, 
                            fitControl, seed){
  ##Combine features and resp
  FeatResp_data <- combine_FeatResp(feats = feats, resp = resp)
  ##Train model
  set.seed(seed)
  mod <- train(response ~ ., data = FeatResp_data,
                method = "ranger",
                metric = "ROC", trControl = fitControl)
  return(mod)
}



##### Make predictions and extract results #####
# Function: Apply models to trainTestSplits to get test-set predictions, accuracy, and final predictions
# Returns: 
## modsPredsTest_df - test-set predictions across trainTestSplits
## AUCdf - AUCs within each test-set
## allPredsDf - all models applied to datFeats, then median of scores is taken
makePreds_fromMods <- function(
  mods,  # List of models, each element corresponds to a train-test split
  trainTestSplits,  # Train-test splits as produced by makeFoldsProtstrat(). Names must match names of mods
  datFeats,   # Dataframe with all features for all observations. Used to make final predictions
  returnAllPredsLst = F # Return all predictions for all models
){
  resultsLst <- list()
  
  ##### Apply to test set #####
  ## Make predictions
  modsPredsTest <- map2(mods, names(mods), function(tempMod, fold){ #Loop over train-test splits
    testDat <- trainTestSplits[[fold]]$testDat
    predictions <- predict(tempMod, newdata = testDat[, -which(colnames(testDat) == "response")], type = "prob")
    predDf <- data.frame("X1" = predictions[, "X1"], "obs" = testDat$response)
    predDf$fold <- fold
    predDf$uniprot_site <- rownames(testDat)
    return(predDf)
  })
  modsPredsTest_df <- purrr::reduce(modsPredsTest, rbind) %>% mutate(Resample = fold)
  resultsLst$modsPredsTest_df <- modsPredsTest_df # Add to resultsLst

  # Calculate AUCs
  AUCs <- modsPredsTest_df %>% split(.$Resample) %>% 
    map(~ROCAUC_fromPredictor(data = ., obsCol = "obs", predCol = "X1")) %>% unlist
  AUCdf <- data.frame("Resample" = names(AUCs), "AUC" = AUCs)
  resultsLst$AUCdf <- AUCdf # Add to resultsLst

  # #Plot ROC
  # ROCdat <- prep_ROC_meansd_fromLst(modsPreds_trainTestforPred)
  # ROCplt <- plotROC_fromData(data = ROCdat, 
  #                  MeanSd = T, linewidth = 0.235*1.5) +
  #                  ggtitle("Repeated test sets")


  ##### Make predictions on all data and combine into score #####
  allPredsLst <- map2(mods, names(mods), function(tempMod, fold){ #Loop over train-test splits
          predictions <- predict(tempMod, newdata = datFeats, type = "prob")
          predDf <- data.frame("X1" = predictions[, "X1"], "uniprot_site" = rownames(datFeats))
          return(predDf)
        })
  allPredsLstDf <- map_dfr(allPredsLst, function(x){x}, .id = "fold")
  #Take median of predictions across sub-sampled models
  allPredsDf <- allPredsLstDf %>% group_by(uniprot_site) %>% summarise(X1 = median(X1, na.rm = T))
  resultsLst$allPredsDf <- allPredsDf # Add to resultsLst
  if (returnAllPredsLst){
    resultsLst$allPredsLst <- allPredsLst
  }
  return(resultsLst)
}

# Function: Apply models to trainTestSplits to get test-set predictions, accuracy, and final predictions
# An example of this in action can be found in 11d_training_pipeline.Rmd
# Returns: 
## modsPredsTest_df - test-set predictions across trainTestSplits
## AUCdf - AUCs within each test-set
## allPredsDf - all models applied to datFeats, then median of scores is taken
## mods - list of models for each train-test split
trainAndMakePreds_logreg <- function(
  trainTestSplits,   # Train-test splits as produced by makeFoldsProtstrat()
  datFeats,  # Dataframe with all features for all observations. Used to make final predictions
  seed = 80,  # Seed for training models
  returnAllPredsLst = F  # Return all predictions for each model?
){
  ##### Train model #####
  mods <- list()
  for (j in 1:length(trainTestSplits)){
    splt <- trainTestSplits[[j]]
    trainDat <- splt$trainDat
    # Train logistic regression
    mod <- modTrain_LogReg_weighted(feats = trainDat[, -which(colnames(trainDat) == "response")], 
                                            resp = trainDat$response, 
                                            fitControl = fitControl_CV, seed = seed + j,
                                            weightFunc = caseWeights_balanced)
      mods[[j]] <- mod
  }
  names(mods) <- names(trainTestSplits)
  ##### Make predictions and get ROC results #####
  resultsLst <- makePreds_fromMods(
    mods = mods,  
    trainTestSplits = trainTestSplits,  
    datFeats = datFeats,
    returnAllPredsLst = returnAllPredsLst
  )
  # Add models to results
  resultsLst$mods <- mods
  return(resultsLst)
}
