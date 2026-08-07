###Last updated:
#20211020

###Background
#Here are all of my functions that perform enrichment
#Fisher's exact
#GSEA
#KSEA

###Packages
library(reshape2)
library(scales)
library(GO.db)
library(rlist)
library(multcomp)
library(ggplot2)
library(tidyverse)
library(purrr)

#####Pathway enrichment#####
####Function that enriches from a list of pathways
#As of 20250116 it works when there are non-unique names in DE_genes and background_genes
#Supply list of relevant pathways, vector of background genes, and vector of DE genes
#Pathway_DE_intersection_threshold: Pathways need this number of DE genes or more to be considered
#intersection_as_unique: if TRUE, only use unique DE_genes for the intersection threshold
pathway_enricher_from_list <- function(background_genes,
                                       DE_genes,
                                       pathways_list,
                                       pathway_DE_intersection_threshold = FALSE,
                                       intersection_as_unique = FALSE,
                                       alternative = "greater"){
  
  ##Make pathways relevant 
  #Trim pathway list so it's only pathways with DE genes
  if (pathway_DE_intersection_threshold == FALSE){
    NULL
  } else {
    relevant_pathways_bool <- NULL
    old_pathways_list <- pathways_list
    for (i in 1:length(old_pathways_list)){
      if (intersection_as_unique){
        relevant_pathways_bool[i] <- length(intersect(DE_genes,
                                                      old_pathways_list[[i]])) >= pathway_DE_intersection_threshold
      } else {
        relevant_pathways_bool[i] <- length(which(DE_genes %in% old_pathways_list[[i]])) >= pathway_DE_intersection_threshold
      }
    }
    
    #Terminate here if no relevant pathways
    if (sum(relevant_pathways_bool) == 0){
      
      return(NULL)
    } else {
      
      pathways_list <- as.list(names(old_pathways_list)[relevant_pathways_bool])
      names(pathways_list) <- names(old_pathways_list)[relevant_pathways_bool]
      for (i in 1:length(pathways_list)){
        
        pathways_list[[i]] <- old_pathways_list[[names(pathways_list)[i]]]
      }
    }
  }
  
  
  
  ##Set up output
  #pval
  #adj_pval
  #Number of DE genes in pathway
  #Number of background genes in pathway
  #Number of total genes in pathway
  output_m <- matrix(NA,
                     nrow = length(pathways_list),
                     ncol = 6)
  rownames(output_m) <- names(pathways_list)
  colnames(output_m) <- c("pval",
                          "adj_pval",
                          "num_DE_genes",
                          "num_DE_genes_in_pathway",
                          "num_background_genes_in_pathway",
                          "num_total_genes_in_pathway")
  
  ##Do stats
  #Do one-sided fisher's exact
  
  
  #Get x, m, n, k
  #m is number of DE genes
  #n is number of background, non-DE genes
  #k is number of background genes in pathway
  #x is number of DE genes in pathway
  
  m <- length(DE_genes)
  n <- length(background_genes) - m
  
  for (i in 1:length(pathways_list)){
    
    temp_pathway <- pathways_list[[i]]
    k <- length(which(background_genes %in% temp_pathway))
    x <- length(which(DE_genes %in% temp_pathway))
    #Run test
    test <- fisher.test(rbind(c(x, k - x),
                              c(m - x, n - k + x)),
                        alternative = alternative)
    
    temp_pval <- test$p.value
    output_m[i, ] <- c(temp_pval,
                       NA,
                       m,
                       x,
                       k,
                       length(temp_pathway))
  }
  
  ##Adjust pvals
  output_m[, 2] <- p.adjust(output_m[, 1],
                            method = "fdr")
  
  ##Sort
  #By adj_pval. Don't order if length == 1
  if (nrow(output_m) > 1){
    output_m <- output_m[order(output_m[, 2],
                               decreasing = FALSE), ]
  }
  
  ##Add enrichment score
  output_m <- as.data.frame(output_m)
  output_m$ES <- (output_m$num_DE_genes_in_pathway/length(DE_genes))/
    (output_m$num_background_genes_in_pathway/length(background_genes))
  output_m$ES_log2 <- log2(output_m$ES)
  
  ##Return
  return(output_m)
}


####Function: Get GO pathways that intersect with a list of genes
#DE_genes: Genes of interest
#gene_label: Nami gnconvention for the genes, which ahs to be read by AnnotationDBI. Default is SYMBOL (e.g. Akt2, Slc2a4)
#organism_database: The AnnotationDBI organism database e.g. org.Mm.eg.db
#ontology: The desired GO ontology e.g. CC
relevant_GOIDS <- function(DE_genes,
                           gene_label = "SYMBOL",
                           organism_database,
                           ontology){
  relevant_GOIDS_df <- AnnotationDbi::select(organism_database,
                                             keys = DE_genes,
                                             columns = "GO",
                                             keytype = gene_label)
  relevant_GOIDS_df <- relevant_GOIDS_df[which(relevant_GOIDS_df$ONTOLOGY == 
                                                 ontology),]
  relevant_GOIDS <- unique(relevant_GOIDS_df$GO)
  
  ###List of all genes for relevant_GOIDS
  relevant_GOIDS_genes_df <- AnnotationDbi::select(organism_database,
                                                   keys = relevant_GOIDS,
                                                   columns = gene_label,
                                                   keytype = "GO")
  ##Make list
  #Make sure no gene name complete duplicates. Can happen because of different evidence levels
  relevant_GOIDS_genes_list <- as.list(unique(relevant_GOIDS_genes_df$GO))
  names(relevant_GOIDS_genes_list) <- unique(relevant_GOIDS_genes_df$GO)
  for (i in 1:length(relevant_GOIDS_genes_list)){
    
    temp_GOID <- names(relevant_GOIDS_genes_list)[i]
    relevant_GOIDS_genes_list[[i]] <- unique(relevant_GOIDS_genes_df[which(relevant_GOIDS_genes_df$GO == temp_GOID), 
                                                                     gene_label])
  }
  return(relevant_GOIDS_genes_list)
}


###Function to perform GO enrichmen given list
#Up to date 20210528
#Supply vector of background genes, vector of DE genes, gene_label (UNIPROT, SYMBOL, etc), organism database (e.g. org.Mm.eg.db), and ontology (e.g. CC)
GO_enricher_perform <- function(background_genes,
                                DE_genes,
                                relevant_GOIDS_genes_list,
                                pathway_DE_intersection_threshold,
                                alternative = "greater",
                                ...){
  ##Run tests
  output_df <- as.data.frame(pathway_enricher_from_list(background_genes = background_genes,
                                                        DE_genes = DE_genes,
                                                        pathways_list = relevant_GOIDS_genes_list,
                                                        pathway_DE_intersection_threshold = pathway_DE_intersection_threshold,
                                                        alternative = alternative,
                                                        ...),
                             stringsasfactors = FALSE)
  
  ##Add GO term and ID
  output_df$GO_term <- AnnotationDbi::select(GO.db,
                                             rownames(output_df),
                                             columns = "TERM",
                                             keytype = "GOID")$TERM
  output_df$GO_id <- rownames(output_df)
  output_df$GO_term <- ifelse(is.na(output_df$GO_term), output_df$GO_id, output_df$GO_term)
  rownames(output_df) <- output_df$GO_term
  
  ##Return
  return(output_df)
}


####Function that enriches for terms of a given GO ontology
#Up to date 20210528
#Supply vector of background genes, vector of DE genes, gene_label (UNIPROT, SYMBOL, etc), organism database (e.g. org.Mm.eg.db), and ontology (e.g. CC)
GO_enricher <- function(background_genes,
                        DE_genes,
                        gene_label = "SYMBOL",
                        organism_database,
                        ontology,
                        pathway_DE_intersection_threshold,
                        alternative = "greater",
                        ...){
  
  ###Get relevant GOIDs
  relevant_GOIDS_genes_list <- relevant_GOIDS(DE_genes = DE_genes,
                                              gene_label = gene_label,
                                              organism_database = organism_database,
                                              ontology = ontology)
  
  
  ##Run tests
  output_df <- GO_enricher_perform(DE_genes = DE_genes,
                                   background_genes = background_genes,
                                   relevant_GOIDS_genes_list = relevant_GOIDS_genes_list,
                                   pathway_DE_intersection_threshold = pathway_DE_intersection_threshold,
                                   alternative = alternative,
                                   ...)
  
  ##Return
  return(output_df)
}



###Fisher's exact test for DE genes and pathway
#Up to date 20210712
#Supply background genes, DE genes, pathway genes, and specifict alternative
fishers_DE_pathway <- function(background_genes,
                               DE_genes,
                               pathway_genes,
                               alternative = "greater"){
  
  ##Run test
  #m: number of DE genes
  #n: number of background, non-DE genes
  #k: number of background genes in pathway
  #x: number of DE genes in pathway
  m <- length(DE_genes)
  n <- length(background_genes) - m
  k <- length(intersect(background_genes,
                        pathway_genes))
  x <- length(intersect(DE_genes,
                        pathway_genes))
  
  test <- fisher.test(rbind(c(x, k - x),
                            c(m - x, n - k + x)),
                      alternative = alternative)
  pval <- test$p.value
  odds_ratio <- test$estimate
  
  ##Fold enrichment
  fold_enr <- (x/m)/(k/length(background_genes))
  
  ##Output
  output <- c(pval,
              x,
              k,
              length(pathway_genes),
              odds_ratio,
              fold_enr)
  output <- as.list(output)
  names(output) <- c("pval",
                     "num_DE_genes_in_pathway",
                     "num_background_genes_in_pathway",
                     "num_pathway_genes",
                     "odds_ratio",
                     "fold_enr")
  return(output)
}