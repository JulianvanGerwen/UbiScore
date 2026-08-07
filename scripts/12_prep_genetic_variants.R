##### Background #####
# Here I prepare genetic variant data e.g. ClinVar 
# Source from the project directory

##### Initialise #####
library(tidyverse)
library(purrr)
load("data/features/ubiFeatProc__1.RData") #Load data so I can filter for sites
home_directory <- ""

##### ClinVar #####
# Load ClinVar data 
clinvarDatRaw <- read_csv(paste0(home_directory, "data/biol_databases/clinvar/clinvarJurgen.csv"))

#Add uniprot_site
clinvarDatRaw$uniprot_site <- strsplit(clinvarDatRaw$variant_id, "/") %>% map(function(vec){
  n <- str_length(vec[2])
  vec[2] <- substr(vec[2], 1, n - 1)
  return(paste(vec, collapse = "_"))
}) %>% unlist

#Subset for uniprot_site in our data
clinvarDatSub <- subset(clinvarDatRaw, uniprot_site %in% ubiFeatProc$uniprot_site)
dim(clinvarDatSub)
length(unique(clinvarDatSub$uniprot_site))
save(clinvarDatSub, file = "output/data/genetics/clinvarDatSub__2.RData")
