##### Background #####
#Here I store scripts for using the ubn omics data


##### Inspecting data #####
# Function: Extract and sort omics FCs for a set of sites
# sites: uniprot_site format
# omicsLst: A list like ubnOmicsComp_hsPRIDE_lst
extractOmicsConds <- function(sites, omicsLst){
    # Replace matrix column names with descriptions
    colnames(omicsLst$matrix) <- map(colnames(omicsLst$matrix), ~subset(omicsLst$meta, expCond_id == .)$combinedDescUID) %>% unlist
    # Fetch and sort FCs
    FClst <- map(setNames(sites, sites), function(site){
        site <- gsub("_[A-Z]", "_", site)
        if (site %in% rownames(omicsLst$matrix)){
            vec <- omicsLst$matrix[site, ]
            vec <- vec[which(!is.na(vec))]
            return(sort(vec, decreasing = T))
        } else {
            return(NULL)
        }  
    })
    return(FClst)
}