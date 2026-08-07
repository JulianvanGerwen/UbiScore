##### Background #####
#Here I store objects for use throughout the project

##### Colours #####
upDownColours <- c("up" = "#d6604d", "down" = "#4393c3")
ubiSiteColour <- "#2172b6"
ubiSiteColour2 <- "#2374b6"
#PSP colours
PSPcolours <- c("notReg" = "#b1b1b1",
                "reg" = "#2071b6",
                "notProtDeg" = "#2071b6",
                "protDeg" = "#e55838",
                "protDegOnly" = "#e55838",
                "hotspot" = "#d3b924")

#PSP colours for outlines. The grey needs to be stronger for this
PSPcoloursLine <- c("notReg" = "#8c8c8c",
                "reg" = "#2071b6",
                "notProtDeg" = "#2071b6",
                "protDeg" = "#e55838",
                "protDegOnly" = "#e55838",
                "hotspot" = "#d3b924")

#Hotspot plot colours
hotspotLineColours <- c("nlogp" = "#417fc1", "enrichment" = "#000000")

#PRIDE colours
prideSetColours <- c("Gold" = "#28265d",
                     "Silver" = "#2c8dcb",
                     "Bronze" = "#88ceeb")

#omics treatments
#omicsTreatColours <- c("Proteasome inhibition" = "#b21f2c",
#                       "DUB inhibition" = "#d87166",
#                       "Neddhilation inhibition" = "#fddac7",
#                       "DNA damage" = "#2167ac",
#                       "Translation inhibition" = "#71afc8",
#                       "ER stress" = "#d1e5ef",
#                       "Protein folding stress" = "#542d88",
#                       "Oxidative stress" = "#867fad",
#                       "Infection" = "#d78ab6",
#                       "Other" = "#ababab")

omicsTreatColours <- c("Proteasome inhibition" = "#9a5aa4",
                       "DUB inhibition" = "#8177ac",
                       "Neddylation inhibition" = "#e8d5e8",
                       "DNA damage" = "#2d8a91",
                       "Translation inhibition" = "#6fc1a1",
                       "ER stress" = "#c8e8e4",
                       "Protein folding stress" = "#4293c3",
                       "Oxidative stress" = "#3a3a3a",
                       "Infection" = "#bababa",
                       "Other" = "#ffffff")



##### Feature names #####
#Map of feature column names to descriptions
featureNames <- list(
    list("PRIDE_gold", "Gold identification confidence", "Proteomics"),
    list("FinalBoxProcNum", "Conservation level", "Evolutionary"),
    list("HotspotNew_bool", "In hotspot", "Evolutionary"),
    list("mutR_am_pathogenicity_all", "AlphaMissense pathogenicity (K-R)", "Other"),
    list("interface_label", "At interface", "Structural"),
    list("NSP_p_q3_H", "Helix probability", "Structural"),
    list("NSP_p_q3_C", "Coil probability", "Structural"),
    list("NSP_rsa", "Relative surface accessibility", "Structural"),
    list("NSP_disorder", "Disorder", "Structural"),
    list("PSP_boolOtherLysMod", "Site contains other PTMs", "Other"),
    list("PSP_numAdj_5", "# nearby PTMs", "Other"),
    list("numQuant_all", "# quantifications", "Proteomics"),
    list("numReg_notPrtsm", "# regulations", "Proteomics"),
    list("UNIPROT_in_region", "In region (uniprot)", "Structural"),
    list("UNIPROT_in_domain", "In domain (uniprot)", "Structural")
)
featureCols <- unlist(map(featureNames, ~.[1]))
featureGroups = c("Evolutionary", "Proteomics", "Structural", "Other")
featureGroupCols <- c("Evolutionary" = "#7d3d89",
                      "Proteomics" = "#20675c",
                      "Structural" = "#42843e",
                      "Other" = "#7f7f7f")
featureGroupFills <- c("Evolutionary" = "#cbbad3",
                      "Proteomics" = "#9dcfc6",
                      "Structural" = "#d3e2ac", 
                      "Other" = "#dddddd")