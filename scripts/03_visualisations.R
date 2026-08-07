#####Background#####
#Here I store code for visualisations


#####Feature exploration#####
#Function: Plot a barplot using summary from catBoolSummariser()
bplotRegSumm <- function(regSummDat, cols, suffix = "enr", xaxis){
  #Prep data
  datLong <- regSummDat %>% pivot_longer(cols = paste(cols, suffix, sep = "_"), names_to = "type", values_to = "perc") %>%
    mutate(type = gsub(paste0("_", suffix, "$"), "", type)) %>% mutate(type = factor(type, levels = cols))
  #Plot
  plt <- ggplot(datLong, aes(x = !!sym(xaxis), y = perc, fill = type, colour = type)) +
    geom_col(position = position_dodge2(), width = 0.8, size = 0.235) +
    comfy_theme(include_xaxis = F) 
  return(plt)
}

#Function: Plot bplotRegSumm for both group percentage and fold enrichment
bplotRegSumm_percAndEnr <- function(regSummDat, cols, xaxis, xlab){
  lst <- list(
    "perc" = bplotRegSumm(regSummDat = regSummDat, cols = cols,  xaxis = xaxis, suffix = "groupPct") +
      labs(x = xlab, y = "% of sites in conservation box"),
    
    "enr" = bplotRegSumm(regSummDat = regSummDat, cols = cols,  xaxis = xaxis, suffix = "enr")+
      geom_hline(yintercept = 1, size = 0.235, linetype = "dashed") +
      labs(x = "Conservation box", y = "Fold enrichment (over all ubi-sites)")
  )
  return(lst)
}


#####Summaries#####
#Function: Barplot from binning
#Group a numerical column into bins and then count number of entries in each bin, in each value of a group column
#The bin breaks denote the starting point of each break. The bins are closest on the left and open on the right. 
#So in the given example, the first two bins are [0, 6), [6, 11)
barFBinAndSummarise <- function(data, bin_breaks = c(0, 6, 11, 21, Inf),
                                groupCol, valCol,
                                breakType = "integer",
                                barWidth = 0.85){
  #Bin data
  summDat <- binAndSummarise(data, 
                             bin_breaks = bin_breaks,
                             groupCol = groupCol, valCol = valCol,
                             breakType = breakType)
  #Plot
  plt <- ggplot(summDat, aes(x = Bin, y = Percentage, fill = group)) +
    geom_col(position = position_dodge(),
             width = barWidth) +
    comfy_theme(include_xaxis = F) +
    labs(y = "% of site group")
  return(plt)
}



#####Omics#####
#Plot correlation between two conditions
corrplotConds <- function(ubnOmicsLst, cond1, cond2, sites_to_label = c()){
  plotDat <- ubnOmicsLst$matrix[, c(cond1, cond2)] %>%
    as.data.frame %>%
    mutate(uniprot_pos = rownames(.))
  plotDat <- plotDat[which(rowSums(is.na(plotDat)) == 0), ]
  r_value <- round(cor(plotDat[[cond1]], plotDat[[cond2]], use = "pairwise.complete.obs"), 2)
  #Set up sites to label
  if (length(sites_to_label) > 0){
    labDat <- subset(plotDat, uniprot_pos %in% sites_to_label)
    labPoint <- geom_point(data = labDat, aes(x = !!sym(cond1), y = !!sym(cond2)), shape = 21, size = 1,
                            stroke = 0.5, colour = "black", fill = NA)
    labText <- geom_text_repel(data = labDat, 
                   aes(x = !!sym(cond1), y = !!sym(cond2), label = uniprot_pos), 
                   size = 1.875,
                   colour = "black")
  } else {
    labPoint <- NULL 
    labText <- NULL
  }
  #Plot
  plt <- ggplot(plotDat, aes(x = !!sym(cond1), y = !!sym(cond2))) +
    geom_hline(yintercept = 0, size = 0.235, colour = "#939393", linetype = "dashed") + 
    geom_vline(xintercept = 0, size = 0.235, colour = "#939393", linetype = "dashed") + 
    stat_smooth(method = "lm", colour = "#cb1a1e", fill = "#cb1a1e", alpha = 0.15, size = 0.235) +
    geom_point(shape = 16, size = 1, alpha = 0.5, colour = "#28275d") + 
    labPoint + labText +
    annotate("text", x = -Inf, y = Inf, label = paste("r =", r_value), 
            hjust = -0.1, vjust = 1.1, size = 1.875*7/5) +
    comfy_theme() +
    labs(x = ubnOmicsLst$meta$combinedDescUID[which(ubnOmicsLst$meta$expCond_id == cond1)], 
        y = ubnOmicsLst$meta$combinedDescUID[which(ubnOmicsLst$meta$expCond_id == cond2)])
  return(plt)
}


##### Experiments #####
#Function: Order an identifier column in long data by hierarchical clustering
order_from_cluster <- function(data, idCol, val = "val", namesFrom){
  #Rename variables
  data$idCol <- data[, idCol]
  data$namesFrom <- data[, namesFrom]
  data$val <- data[, val]
  ##Pivot to wide
  dataWide <- pivot_wider(data[, c("idCol", "namesFrom", "val")], names_from = namesFrom, values_from = val) %>%
    as.data.frame
  rownames(dataWide) <- dataWide[, 1]
  dataWide <- dataWide[, -1]
  ##Cluster
  dist_m <- dist(dataWide)
  hclust <- hclust(dist_m)
  dataWide <- dataWide[hclust$order, ]
  return(rownames(dataWide))
}


#heatmap of sscores
hmap_sscores <- function(sscoreDat, sigThresh = 0.01, tileBorderSize = 0,
                         xaxis = "replicateID", yaxis = "condition", clusterY = T){
  ## Cluster y-axis
  if (clusterY & length(unique(sscoreDat[, xaxis])) > 1){
    yClust <- order_from_cluster(data = sscoreDat, idCol = yaxis, val = "score", namesFrom = xaxis)
    sscoreDat[, yaxis] <- factor(sscoreDat[, yaxis], levels = rev(yClust))
    }
  ##Set up colouring
  maxVal <- max(abs(sscoreDat$score), na.rm = TRUE)
  gradient <- scale_fill_gradientn(colors = rev(brewer.pal(9, "RdBu")),
                        limits = c(-maxVal, maxVal),
                        name = "s-score")
  ##Set up dots
  sscoreDat$sig <- sscoreDat$qvalue < sigThresh
  pointGeom <- geom_point(data = subset(sscoreDat, sig == T),
                          aes(x = !!sym(xaxis), y = !!sym(yaxis)),
                          size = 0.5, shape = 16, colour = "black")
  ##Make plot
  plt <- ggplot(sscoreDat, aes(x = !!sym(xaxis), y = !!sym(yaxis), fill = score)) +
    geom_tile(colour = "#a5a5a5",
              size = tileBorderSize) + 
    pointGeom +
    comfy_theme(include_xaxis = F, include_yaxis = F,
                rotate_x_text = T) +
    theme(axis.text.y = element_text(colour = "black", size = 6, hjust = 1)) +
    coord_fixed(ratio = 1/2) +
    gradient
  return(plt)
}



#General heatmap
#Function: Create a heatmap from numerical data with x and y axes clustered
#vars: A named list of variables following the given template
#colourScale: The input to "values" in "scale_fill_gradientn"
#tileBorderSize: Used for outlining cells. 0.05 is a good value
hmapClusteredGeneral <- function(data,
                                 vars = list(
                                   "x" = "x",
                                   "y" = "y",
                                   "val" = "val",
                                   "pval" = "pval"
                                 ),
                                 clusterY = T,
                                 clusterX = T,
                                 colours = c("#ffff03", "#0d0d0d", "#18b7e8"),
                                 colourScale = c(-1, -0.33, 0, 0.33, 1),
                                 includeXText = FALSE,
                                 includeYText = FALSE,
                                 tileBorderSize = FALSE,
                                 pval_cutoff = FALSE){
  data <- as.data.frame(data)
  #Rename variables
  for (i in 1:length(vars)){
    data[, names(vars)[i]] <- data[, vars[[i]]]
  }
  #Cluster variables
  if (clusterY){
    yClust <- order_from_cluster(data = data, idCol = "y", namesFrom = "x")
    data <- mutate(data, y = factor(y, levels = rev(yClust)))
    }
  if (clusterX){
    xClust <- order_from_cluster(data = data, idCol = "x", namesFrom = "y")
    data <- mutate(data, x = factor(x, levels = rev(xClust)))
    }
  #Significance dots
  if (pval_cutoff == FALSE){
    point <- NULL
  } else {
    data$sig <- data$pval < pval_cutoff
    point <- geom_point(data = data[data$sig == TRUE, ],
                        aes(x = x, y = y),
                        shape = 16, size = 0.5, colour = "red")
  }
  #Colour scale
  maxVal <- max(abs(data$val)) #Get max for limits
  #Theme
  if (!includeXText){
    xtheme <- theme(axis.text.x = element_blank())
  } else {
    xtheme <- NULL
  }
  if (!includeYText){
    ytheme <- theme(axis.text.y = element_blank())
  } else {
    ytheme <- NULL
  }
  #Tile geom
  if (!tileBorderSize){
    tile <- geom_tile()
  } else {
    tile <- geom_tile(colour = "#a5a5a5", size = tileBorderSize)
  }
  #Plot
  outputPlot <- ggplot(data, aes(x = x, y = y, fill = val)) +
    tile + point +
    comfy_theme(rotate_x_text = TRUE, x_text_angle = 90,
                include_xaxis = F, include_yaxis = F) +
    xtheme + ytheme +
    labs(x = vars$x, y = vars$y) +
    scale_fill_gradientn(colors = colours, limits = c(-maxVal, maxVal),
                         values = rescale(colourScale),
                         name = vars$val)
  return(outputPlot)
}

##### Functional score #####
# Function: Plot a lollypop plot of a protein, including domains
# Example usage:
# plotDat <- subset(ubiFeatPreds, gene == "ELAVL1") %>% 
#     mutate(label = pos == 320)
# plot_lollipop(data = plotDat, 
#     domains = list(
#     "RRM1" = c(20, 98),
#     "RRM2" = c(106, 186),
#     "RRM3" = c(244, 322)
#     ), 
#     aes_colour = "label", aes_size = "label",
#     colour_values = setNames(c(ubiSiteColour, "#69a4c4"), c(T, F)), size_values = setNames(c(2.5, 1.5), c(T, F)))
plot_lollipop <- function(data, domains, colour_values, size_values, aes_colour, aes_size, ycol = "functional_score",
                          sitesToLabel = c()) {
    # Get protein length
    protLength <- data$proteinLength[1]

    # Convert domains list to a data frame for plotting rectangles
    domain_df <- tibble::tibble(
        domain = names(domains),
        xmin = sapply(domains, function(x) x[1]),
        xmax = sapply(domains, function(x) x[2])
    )

    # Parameters
    rectHeight <- 0.05 * 1.25 * 1.25

    # Set up label
    if (length(sitesToLabel) > 0){
      data$toLabel <- data$uniprot_site %in% sitesToLabel
      labelPlt <- geom_text(
        data = subset(data, toLabel == T), aes(x = pos, y = !!sym(ycol) + 0.2, colour = !!sym(aes_colour), label = site),
        size = 7/5*1.875, inherit.aes = F)
    } else {
      labelPlt <- NULL
    }

    # Make lollipop plot
    plt <- ggplot(data, aes(x = pos, y = !!sym(ycol), colour = !!sym(aes_colour), size = !!sym(aes_size))) +
        # Background lines
        geom_hline(yintercept = 0, size = 0.235, alpha = 1) +
        geom_hline(yintercept = c(0.5, 1), size = 0.235, alpha = 0.3) +
        # Lollipop
        geom_segment(aes(x = pos, xend = pos, y = 0, yend = !!sym(ycol)), size = 0.235, colour = "#353535") +
        geom_point(shape = 16) + labelPlt +
        # Background rectangle
        geom_rect(aes(xmin = 0, xmax = protLength, ymin = -rectHeight, ymax = rectHeight),
                  fill = "#e9e9ea", colour = "#bbbbba", linewidth = 0.235, inherit.aes = FALSE) +
        # Add rectangles for domains centered around y = 0
        geom_rect(data = domain_df, 
                  aes(xmin = xmin, xmax = xmax, ymin = -rectHeight, ymax = rectHeight),
                  fill = "#cfcfd1", colour = "#bbbbba", linewidth = 0.235, inherit.aes = FALSE) +
        # Plot text inside the rectangles with the domain name
        geom_text(data = domain_df, 
                  aes(x = (xmin + xmax) / 2, y = 0, label = domain),
                  inherit.aes = FALSE, size = 1.875 * 5 / 5, vjust = 0.5, hjust = 0.5, colour = "#555655") +
        comfy_theme(include_xaxis = TRUE, include_yaxis = FALSE, include_legend = FALSE) + 
        theme(axis.line.x = element_blank()) +
        # Scale
        scale_colour_manual(values = colour_values) +
        scale_size_manual(values = size_values) +
        # More elements
        xlim(c(0, protLength)) +
        scale_y_continuous(breaks = seq(0, 1, by = 0.5)) +
        labs(x = "Protein sequence", y = "Ubi-site positional importance score") +
        ggtitle(data$gene[1])
    return(plt)
}


# Function: Plot barplot of coefficients with standard errors, coloured and ordered
# coeffs: List of coefficient values coming from a list of multiple (resampled) models
bplot_coeffs <- function(coeffs){
    # Turn into data frame
    coeffNames <- names(coeffs[[1]])
    modelCoeffsRaw <- map_dfr(coeffs, function(coeffs){
            data.frame("coeff" = coeffNames, "val" = coeffs[coeffNames])
    }, .id = "rep")
    # Drop intercept and flip sign
    modelCoeffs <- subset(modelCoeffsRaw, coeff != "(Intercept)") %>% 
        mutate(val = -val)
        
    # Summarise mean and standard deviation of "val" for each "coeff"
    modelCoeffs_summary <- modelCoeffs %>%
        group_by(coeff) %>%
        summarise(mean_val = mean(val, na.rm = TRUE), sd_val = sd(val, na.rm = TRUE)) %>%
        arrange(mean_val) %>% 
        mutate(coeff = factor(coeff, levels = unique(coeff))) %>% 
        mutate(featureGroup = map(coeff, function(lev){
            bool <- map(featureNames, ~.[[1]] == lev) %>% unlist
            ind <- which(bool == TRUE)
            if (length(ind) > 0){
                return(featureNames[[ind[1]]][[3]])
            } else {
                return(NA)
            }
        }) %>% unlist) %>% 
        mutate(featureGroup = factor(featureGroup, levels = featureGroups))
    levels(modelCoeffs_summary$coeff) <- map(levels(modelCoeffs_summary$coeff), function(lev){
        bool <- map(featureNames, ~.[[1]] == lev) %>% unlist
        return(featureNames[bool][[1]][[2]])
    }) %>% unlist 

    # Barplot, order and colour by feature groups
    plt <- modelCoeffs_summary %>%
        arrange(desc(featureGroup), mean_val) %>%
        mutate(coeff = factor(coeff, levels = unique(coeff))) %>%
    ggplot(., aes(x = mean_val, y = coeff, colour = featureGroup, fill = featureGroup)) +
        geom_bar(stat = "identity", width = 0.8, linewidth = 0.235) +
        geom_errorbar(aes(xmin = mean_val - sd_val, xmax = mean_val + sd_val), 
                                    width = 0.2, color = "black", size = 0.235) +
        scale_colour_manual(values = featureGroupCols) +
        scale_fill_manual(values = featureGroupFills) +
        labs(x = "Feature weight") +
        comfy_theme(include_legend = F) + 
        theme(axis.title.y = element_blank(), axis.text.y = element_text(colour = "black", size = 6, hjust = 1))
    return(plt)
}


# Function: Barplot of mean and sd for some metric across model features
# Will rename features using nice names and colour them using feature groups
# featCol: The column that indicates the features
bar_meanSd_features <- function(
    data, 
    meanCol = "mean", sdCol = "sd",
    featCol = "model", orderByMean = T,
    xLabel = "AUC", includeVline = T){
    # Assign columns
    data$mean <- data[[meanCol]]
    data$sd <- data[[sdCol]]
    data$feature <- data[[featCol]]
    # Assign feature groups and names
    data <- data %>%
            mutate(feature = factor(feature, levels = unique(feature))) %>% 
            # Get feature groups
            mutate(featureGroup = map(feature, function(lev){
                bool <- map(featureNames, ~.[[1]] == lev) %>% unlist
                ind <- which(bool == TRUE)
                if (length(ind) > 0){
                    return(featureNames[[ind[1]]][[3]])
                } else {
                    return(NA)
                }
            }) %>% unlist) %>% 
            mutate(featureGroup = factor(featureGroup, levels = featureGroups))
    # Reassign levels using proper feature names
    levels(data$feature) <- map(levels(data$feature), function(lev){
        bool <- map(featureNames, ~.[[1]] == lev) %>% unlist
        ind <- which(bool == TRUE)
        if (length(ind) > 0){
            return(featureNames[[ind[1]]][[2]])
        } else {
            return(lev)
        }
    }) %>% unlist 
    # Order if desired
    if (orderByMean){
        data <- data %>% .[order(.$mean, decreasing = F), ] %>% 
        mutate(feature = factor(feature, levels = unique(as.character(feature))))
    }
    # Add vertical lines for AUC if desired
    if (includeVline){
      vlinePlot <- geom_vline(xintercept = c(0.5, 0.6, 0.7), linewidth = 0.235, alpha = 0.5) 
    } else {vlinePlot <- NULL}
    # Make plot
    plt <- ggplot(data, aes(x = mean, y = feature, colour = featureGroup, fill = featureGroup)) +
            vlinePlot +
            geom_bar(stat = "identity", width = 0.8, linewidth = 0.235) +
            geom_errorbar(aes(xmin = mean - sd, xmax = mean + sd), 
                                        width = 0.2, color = "black", size = 0.235) +
            scale_colour_manual(values = featureGroupCols) +
            scale_fill_manual(values = featureGroupFills) +
            comfy_theme(include_legend = F) + 
            labs(x = xLabel) +
            theme(axis.text.y = element_text(colour = "black", size = 7, hjust = 1)) 
    return(plt)
}
