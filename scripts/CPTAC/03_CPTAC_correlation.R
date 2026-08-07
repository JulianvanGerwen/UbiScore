###### Background #####
# Functions to correlate CPTAC ubiquitylome sites against kinase/TF activity,
# regressing out protein abundance, then annotate and filter the resulting hits.
# Source this from the project directory, not from its own directory

##### Initialise #####
library(tidyverse)
library(purrr)

##### Functions #####

#Regress each column in `cols` against `prot`, appending "<col>_resid" columns
regress_out_protein <- function(siteDat, cols = c("site", "activity")){
    for (col in cols){
        siteDat$y <- siteDat[, col]
        resCol <- paste(col, "resid", sep = "_")
        siteDat[, resCol] <- NA
        if (sum(complete.cases(siteDat[, c("y", "prot")])) > 2){
            lm_site <- lm(y ~ prot, data = siteDat)
            siteDat[as.numeric(names(residuals(lm_site))), resCol] <- residuals(lm_site)
        }
        siteDat <- dplyr::select(siteDat, -c("y"))
    }
    return(siteDat)
}

#Correlate two vectors on their mutual complete observations
correlate_vectors <- function(corVec1, corVec2){
    common_idx <- complete.cases(corVec1, corVec2)
    n_obs <- sum(common_idx)
    if (n_obs > 2){
        cor_test <- cor.test(corVec1[common_idx], corVec2[common_idx])
        r_val <- as.numeric(cor_test$estimate)
        p_val <- cor_test$p.value
    } else {
        r_val <- NA
        p_val <- NA
    }
    return(c("r" = r_val, "p" = p_val, "n_obs" = n_obs))
}

#Correlate a single activity matrix against ubi-site quantification for one CPTAC dataset
run_activity_correlation_single <- function(activMat, CPTAC_ubi_dat, CPTAC_ubi_numColsTum, ubiFeatPreds_byGene, CPTACprot){
    regrDatLst_sub <- list()
    #Only correlate over samples quantified in the ubi-data, the activity matrix and the proteomics
    sampleCols <- intersect(CPTAC_ubi_numColsTum, intersect(colnames(activMat), colnames(CPTACprot)))
    activSites <- subset(ubiFeatPreds_byGene, gene %in% rownames(activMat)) %>% #Subset for our ubi-sites
                    subset(gene_site %in% CPTAC_ubi_dat$gene_site)
    result <- data.frame(gene = character(),
                        gene_site = character(),
                        r_value = numeric(),
                        p_value = numeric(),
                        p_value_adj = numeric(),
                        n_obs = integer(),
                        r_value_resid = numeric(),
                        p_value_resid = numeric(),
                        p_value_adj_resid = numeric(),
                        n_obs_resid = integer(),
                        stringsAsFactors = FALSE)

    for (i in 1:nrow(activSites)){
        gene <- activSites$gene[i]
        gene_site <- activSites$gene_site[i]

        # Extract vec1
        vec1 <- CPTAC_ubi_dat[which(CPTAC_ubi_dat$gene_site == gene_site), sampleCols]
        if(nrow(vec1) > 1){
            vec1 <- as.numeric(apply(vec1, 2, mean, na.rm = TRUE))
        } else {
            vec1 <- as.numeric(vec1)
        }
        # Extract vec2
        vec2 <- as.numeric(activMat[gene, sampleCols])
        # Extract protein abundance
        vec3 <- as.numeric(CPTACprot[gene, sampleCols])
        # Assemble data and regress out protein abundance
        siteDat <- data.frame("site" = vec1, "activity" = vec2, "prot" = vec3)
        siteDat <- regress_out_protein(siteDat, cols = c("site", "activity"))
        regrDatLst_sub[[i]] <- siteDat

        #Correlate raw values, and residuals after protein normalisation
        resLst <- list(
            "raw" = correlate_vectors(siteDat$site, siteDat$activity),
            "resid" = correlate_vectors(siteDat$site_resid, siteDat$activity_resid)
        )

        # Append the results
        result <- rbind(result, data.frame(gene = gene,
                                            gene_site = gene_site,
                                            r_value = resLst$raw["r"],
                                            p_value = resLst$raw["p"],
                                            p_value_adj = NA,
                                            n_obs = resLst$raw["n_obs"],
                                            r_value_resid = resLst$resid["r"],
                                            p_value_resid = resLst$resid["p"],
                                            p_value_adj_resid = NA,
                                            n_obs_resid = resLst$resid["n_obs"],
                                            stringsAsFactors = FALSE))
    }
    names(regrDatLst_sub) <- activSites$gene_site
    return(list("regrDat" = regrDatLst_sub, "corrRes" = result))
}

#Run the correlation across every activity matrix in activMatLst for one CPTAC dataset
run_activity_correlation_all <- function(activMatLst, CPTAC_ubi_dat, CPTAC_ubi_numColsTum, ubiFeatPreds_byGene, CPTACprot){
    regrDatLst <- list()
    activCorrResLst <- list()
    for (type in names(activMatLst)){
        res <- run_activity_correlation_single(
            activMat = activMatLst[[type]],
            CPTAC_ubi_dat = CPTAC_ubi_dat,
            CPTAC_ubi_numColsTum = CPTAC_ubi_numColsTum,
            ubiFeatPreds_byGene = ubiFeatPreds_byGene,
            CPTACprot = CPTACprot
        )
        regrDatLst[[type]] <- res$regrDat
        activCorrResLst[[type]] <- res$corrRes
    }
    return(list("regrDatLst" = regrDatLst, "activCorrResLst" = activCorrResLst))
}

#Map in functional score, PSP and hotspot annotation, and compute (adjusted) p-values
annotate_correlation_hits <- function(activCorrResLst, ubiFeatPreds_byGene, hotspotSites){
    forJoin <- left_join(ubiFeatPreds_byGene, hotspotSites[, c("uniprot_site", "hotspot_id")], by = "uniprot_site")
    activCorrResLst <- map(activCorrResLst, ~left_join(., forJoin[, c("gene_site", "functional_score",
                                                        "PSP_reg", colnames(forJoin)[grep("^PSPreg", colnames(forJoin))],
                                                        "HotspotNew_bool", "hotspot_id", "BortMG_meanlog2FC")],
                                                        by = "gene_site"))
    activCorrResLst <- map(activCorrResLst, function(dat){
        dat$nlogp <- -log10(dat$p_value)
        dat$p_value_adj <- p.adjust(dat$p_value, method = "fdr")
        dat$p_value_adj_resid <- p.adjust(dat$p_value_resid, method = "fdr")
        return(dat)
    })
    return(activCorrResLst)
}

#Filter annotated correlation results for significant, biologically-consistent hits
filter_correlation_hits <- function(activCorrResLst){
    activCorrResLst_filt <- map(activCorrResLst, function(dat){
        ##Initial filter
        dat <- dat %>%
            mutate(noPrtsmUpreg = is.na(BortMG_meanlog2FC) | BortMG_meanlog2FC < 0) %>%
            dplyr::relocate(c("hotspot_id", "BortMG_meanlog2FC", "functional_score"))  %>%
            subset(p_value_adj < 0.05 & noPrtsmUpreg == T) %>%
            .[order(.$functional_score, decreasing = T), ]
        ##If positively correlated, only retain when regressed p < 0.05
        dat <- mutate(dat, posCorRetain = r_value < 0 | p_value_resid < 0.05) %>%
                subset(posCorRetain == T)
        return(dat)
    })
    return(activCorrResLst_filt)
}

#Full pipeline for one CPTAC dataset: correlate against every activity matrix, annotate, then filter for hits
run_CPTAC_correlation_pipeline <- function(CPTAC_ubi_dat, CPTAC_ubi_numColsTum, activMatLst, ubiFeatPreds_byGene, CPTACprot, hotspotSites){
    corrRes <- run_activity_correlation_all(
        activMatLst = activMatLst,
        CPTAC_ubi_dat = CPTAC_ubi_dat,
        CPTAC_ubi_numColsTum = CPTAC_ubi_numColsTum,
        ubiFeatPreds_byGene = ubiFeatPreds_byGene,
        CPTACprot = CPTACprot
    )
    raw <- annotate_correlation_hits(corrRes$activCorrResLst, ubiFeatPreds_byGene, hotspotSites)
    filtered <- filter_correlation_hits(raw)
    return(list("regrDatLst" = corrRes$regrDatLst, "raw" = raw, "filtered" = filtered))
}
