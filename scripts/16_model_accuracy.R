##### Background #####
# This script contains auxiliary functions that I use to assess the accuracy of our predictions


##### Functions #####
# Function: Calculates AUC within proteins for data with "obs", "Resample", and a column to be used as a numerical predictor
calculateAUC_perProtein <- function(dat, predCol = "X1"){
	# Prepare data
	dat <- as.data.frame(dat)
	dat$predCol <- dat[, predCol]
	# Filter for quantification and both positive and negative cases
	obsFrac <- dat %>% subset(!is.na(predCol)) %>%
		group_by(uniprot) %>% summarise(frac = mean(obs == "X1"))
	unips_wBoth <- subset(obsFrac, frac > 0 & frac < 1)$uniprot # Get only proteins with positive and negative cases
	dat <- subset(dat, uniprot %in% unips_wBoth)
	AUC_perProt <- dat %>% 
		group_by(Resample, uniprot) %>%
			summarise(
				AUC = ROCAUC_fromVecs(obs, predCol),
				num_X1 = sum(obs == "X1"),
				num_X0 = sum(obs == "X0"),
				.groups = "drop"
			)
	return(AUC_perProt)
}

# Function: For data with "obs", "Resample", and a column to be used as a numerical predictor
# Calculates AUC across all sites, within proteins, and across all sites using the protein-level average of predCol
calculateAUC_withinAndBetweenProts <- function(dat, predCol = "X1"){
	## AUC overall
	AUCoverall <- dat %>% group_by(Resample) %>% 
		summarise(AUC = ROCAUC_fromVecs(obs, !!sym(predCol)), .groups = "drop", 
					numX1 = sum(obs == "X1"), numX0 = sum(obs == "X0"))
	## AUC within proteins
	AUCperProt <- calculateAUC_perProtein(dat, predCol = predCol)
	## AUC using protein-level mean
	protLevelMean <- dat %>% group_by(Resample, uniprot) %>% 
		summarise(featMean = mean(!!sym(predCol)))
	protLevelMean <- full_join(dat[, c("obs", "uniprot", "uniprot_site", "Resample")], protLevelMean, by = c("uniprot", "Resample"))
	AUCprotMean <- protLevelMean %>% group_by(Resample) %>% 
		summarise(AUC = ROCAUC_fromVecs(obs, featMean), .groups = "drop", 
					numX1 = sum(obs == "X1"), numX0 = sum(obs == "X0"))
	## Put into list
	resLst <- list("AUCoverall" = AUCoverall, "AUCperProt" = AUCperProt, "AUCprotMean" = AUCprotMean)
	return(resLst)
}

# Function: For data with "Resample" and a column to be used as a numerical predictor
# Calculates MAD across all sites, within proteins, and across all sites using the protein-level average of predCol
calculateMAD_withinAndBetweenProts <- function(dat, predCol = "X1", removeSingles = T){
	dat <- as.data.frame(dat)
	## MAD overall
	MADoverall <- dat %>% group_by(Resample) %>% 
		summarise(MAD = mad(!!sym(predCol), na.rm = T), .groups = "drop", 
					numX1 = sum(obs == "X1"), numX0 = sum(obs == "X0"))
	## MAD within proteins
	MADperProt <- dat %>% group_by(Resample, uniprot) %>% 
	summarise(MAD = mad(!!sym(predCol)), num_X1 = sum(obs == "X1"), num_X0 = sum(obs == "X0"), numTotal = n(), .groups = "drop")
	### Remove proteins wiht only one site if desired
	if (removeSingles){
		MADperProt <- subset(MADperProt, numTotal > 1)
	}
	## MAD using protein-level mean
	protLevelMean <- dat %>% group_by(Resample, uniprot) %>% 
		summarise(featMean = mean(!!sym(predCol)))
	protLevelMean <- full_join(dat[, c("obs", "uniprot", "uniprot_site", "Resample")], protLevelMean, by = c("uniprot", "Resample"))
	MADprotMean <- protLevelMean %>% group_by(Resample) %>% 
		summarise(MAD = mad(featMean, na.rm = T), .groups = "drop", 
					numX1 = sum(obs == "X1"), numX0 = sum(obs == "X0"))
	## Put into list
	resLst <- list("MADoverall" = MADoverall, "MADperProt" = MADperProt, "MADprotMean" = MADprotMean)
	return(resLst)
}


# Function: Combine results from calculateMAD_withinAndBetweenProts or calculateMAD_withinAndBetweenProts,
# by taking the mean within proteins
summarise_AUC_withinAndBetweenProts <- function(resLst, metric = "AUC"){
	lstName <- paste0(metric, "perProt")
	resLst[[lstName]] <- resLst[[lstName]] %>% group_by(Resample) %>% 
		summarise(AUC = mean(!!sym(metric)), num_X1 = mean(num_X1), num_X0 = mean(num_X0), .groups = "drop")
	colnames(resLst[[lstName]])[which(colnames(resLst[[lstName]]) == "AUC")] <- metric 
	resDf <- map_dfr(resLst, function(x){x[, c("Resample", metric)]}, .id = "Regime")
	return(resDf)
}

# Function: Calculate FPR, FNR, Precision, and F1 score using caret::confusionMatrix
# Both obs and pred should be factors with levels "X1" (positive) and "X0" (negative)
calculateClassMetrics <- function(dat, obsCol = "obs", predCol = "pred"){
	obs <- factor(dat[[obsCol]], levels = c("X1", "X0"))
	pred <- factor(dat[[predCol]], levels = c("X1", "X0"))
	
	cm <- confusionMatrix(pred, obs, positive = "X1")
	
	return(data.frame(
		TP = cm$table["X1", "X1"],
		TN = cm$table["X0", "X0"],
		FP = cm$table["X1", "X0"],
		FN = cm$table["X0", "X1"],
		Sensitivity = cm$byClass["Sensitivity"],  # Recall / TPR
		Specificity = cm$byClass["Specificity"],  # TNR = 1 - FPR
		Precision = cm$byClass["Precision"],
		F1 = cm$byClass["F1"],
		FPR = 1 - cm$byClass["Specificity"],
		FNR = 1 - cm$byClass["Sensitivity"]
	))
}