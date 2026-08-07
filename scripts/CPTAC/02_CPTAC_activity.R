###### Background #####
# Here I estimate kinase and transcription factor activity in CPTAC data
# Source this from the project directory, not from its own directory

##### Initialise #####
library(tidyverse)
library(purrr)
library(dorothea)
library(viper)

home_directory <- ""


##### Kinase activity estimation#####
#Set parameters
minSubNum <- 5
invivoOnly <- T

#Functions
runZtest_matrix <- function(
                      phosMat, #Matrix where columns are conditions
                      phosNames, #rownames of phosMat
                      kinSubLst, #A list of substrates of the kinase of interest
                      min_sub_num = 5 #Number of quantified substrates for analysis to run
){
  #Set up
  rownames(phosMat) <- phosNames
  z_m <- p_m <- activ_m <- matrix(NA, ncol = ncol(phosMat), nrow = length(kinSubLst),
                                dimnames = list(names(kinSubLst), colnames(phosMat)))
  #Calculate Global stats
  muT <- apply(phosMat, 2, FUN = "mean", na.rm = T)
  sdT <- apply(phosMat, 2, FUN = "sd", na.rm = T)
  for (kinase in names(kinSubLst)){
    kinSubs <- kinSubLst[[kinase]]
    phosMat_KinSub <- phosMat[intersect(rownames(phosMat), kinSubs), ]
    Sm <- apply(phosMat_KinSub, 2, FUN = "mean", na.rm = T)
    m <- apply(phosMat_KinSub, 2, FUN = function(x){sum(!is.na(x))})
    ##z-scores
    z_score <- ((Sm - muT) * sqrt(m))/sdT
    ##Remove those without enough substrates
    z_score[which(m < min_sub_num)] <- NA
    p_value <- 2*pnorm(-abs(z_score))
    sign <- sign(z_score)
    activity <- -log10(p_value)*sign
    ##Assign to matrices
    z_m[kinase, ] <- z_score
    p_m[kinase, ] <- p_value
    activ_m[kinase, ] <- activity
  }
  return(list("z" = z_m, "p" = p_m, "activity" = activ_m))
}


#Prep CPTAC data
CPTACphos <- read_csv("output/data/CPTAC/phosphoproteomics_merged.csv") %>% as.data.frame
rownames(CPTACphos) <- CPTACphos$gene_site
CPTACphos <- CPTACphos[, -1]

#Prep PhosphositePlus data
PSP_kin_raw <- read_delim(paste(home_directory, "data/biol_databases/phosphosite_plus/DATE/Kinase_Substrate_Dataset", sep = "")) %>%
  as.data.frame
PSP_kin_proc <- PSP_kin_raw %>%
  dplyr::rename("site" = "SUB_MOD_RSD") %>%
  mutate(AA = strsplit(site, "") %>% map(~.[1]) %>% unlist) %>%
  mutate(uniprot_site = paste(SUB_ACC_ID, site, sep = "_")) %>%
  subset(AA %in% c("S", "T", "Y")) %>%
  mutate(gene_site = paste(SUB_GENE, site, sep = "_")) %>%
  subset(SUB_ORGANISM == "human") %>%
  mutate(GENE = toupper(GENE)) #Use kinase gene symbols


#Make list of kinase substrates
kinSubDat <- PSP_kin_proc 
if (invivoOnly){
  kinSubDat <- subset(kinSubDat, IN_VIVO_RXN == "X")
}
kinSubLst <- kinSubDat %>% split(.$GENE) %>% map(~unique(.$gene_site))
##Subset for sites in CPTAC
kinSubLst <- map(kinSubLst, ~intersect(., rownames(CPTACphos)))
kinSubQuants <- map(kinSubLst, length) %>% unlist
sort(kinSubQuants, decreasing = T)
##Filter for kinases with >= k substrates
filtKins <- names(kinSubQuants)[which(kinSubQuants >= minSubNum)]
kinSubLst <- kinSubLst[filtKins]

#Run activity estimation
activ_res <- runZtest_matrix(
                      phosMat = CPTACphos, #Matrix where columns are conditions
                      phosNames = rownames(CPTACphos), #rownames of phosMat
                      kinSubLst = kinSubLst, #A list of substrates of the kinase of interest
                      min_sub_num = minSubNum #Number of quantified substrates for analysis to run
)
z_m <- activ_res$z
p_m <- activ_res$p
activ_m <- activ_res$activity
write.csv(activ_m, file = "output/data/CPTAC/kinase_activity_m.csv")

##### Kinase activity estimation with VIPER #####
#Convert kinase substrates to dorothea format
kinSubsDoro <- map2(kinSubLst, names(kinSubLst), function(vec, kinase){
  dat <- data.frame(
    tf = kinase,
    confidence = "A",
    target = vec,
    mor = rep(1, length(vec))
  )
  return(dat)
}) %>% purrr::reduce(rbind)
dorothea_viper <- df2regulon(kinSubsDoro[,c(3,1,4)]) #Convert to viper format
#Run viper
kinaseActivityViper <- viper(eset = CPTACphos, regulon = dorothea_viper,
                                  minsize = minSubNum, adaptive.size = F, eset.filter = F)
write.csv(kinaseActivityViper, file = "output/data/CPTAC/kinase_activity_VIPER_m.csv")



##### Transcription factor activity #####
#Get all TFs with ubi sites
TFs_all <- subset(dorothea_hs, confidence %in% c("A","B","C"))$tf %>% unique
ubiSitesGenes <- read_csv("output/data/ubiSitesGenes.csv") %>% .[, -1]
TFs_wUbi <- intersect(TFs_all, ubiSitesGenes$gene)

#Load CPTAC data
CPTACtrans <- read_csv("output/data/CPTAC/transcriptomics_merged.csv") %>% as.data.frame
CPTACtransIdentif <- CPTACtrans[, c("Name", "Database_ID")] %>% dplyr::rename("gene" = "Name", "ENSG" = "Database_ID")
rownames(CPTACtrans) <- CPTACtrans$Name
CPTACtrans <- CPTACtrans[, -c(1, 2)]

#Predict activity
dorothea <- as.data.frame(dorothea_hs[dorothea_hs$confidence %in% c("A","B","C"),]) %>% #Get TF target info
  subset(tf %in% TFs_wUbi)
dorothea_viper <- df2regulon(dorothea[,c(3,1,4)]) #Convert to viper format
TFactivity <- as.data.frame(viper(eset = CPTACtrans, regulon = dorothea_viper, #Run TF activity inference 
                                  minsize = 5, adaptive.size = F, eset.filter = F))
write.csv(TFactivity, file = "output/data/CPTAC/TFactivity_m.csv")