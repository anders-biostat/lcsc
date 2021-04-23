
###############################################################################

# Load library

library(sparseMatrixStats)
library(Matrix)

###############################################################################

# Load Data

load_flag = TRUE

if (load_flag) {
  # windows
  dir <- "~/3_r_projects/Anders_Lab_Local/22-04-21_processing/data/"
  # ubuntu
  #dir <- "~/data/own_processed"
  load(paste(dir,"/counts.RData", sep=""))
  load(paste(dir,"/fracs.RData", sep=""))
  load(paste(dir,"/cells.RData", sep=""))
  load(paste(dir,"/sample_info.RData", sep=""))
  load(paste(dir,"/pcx.RData", sep=""))
  load(paste(dir,"/pcs.RData", sep=""))
  load(paste(dir,"/tsne_pcx.RData", sep=""))
  load(paste(dir,"/tsne_pcs.RData", sep=""))
  load(paste(dir,"/umap_pcx.RData", sep=""))
  load(paste(dir,"/umap_pcs.RData", sep=""))
}

# Not implemented yet to compare PCS and PCX in the program dynamically
tsne <- tsne_pcx
umap <- umap_pcx


###############################################################################

if (load_flag) {
  library(tidyverse)
    # Get the reference data
  ref <- read.table("./3_covid_balf/all.cell.annotation.meta.txt", header=TRUE) %>% 
    mutate(colname = str_c(sample_new, "-", str_replace(ID, "_\\d*", ""))) %>%
    rename(sample_old = sample)
  
  # Join the ref with our cell data frame
  cells %>% left_join(ref, by="colname") -> cells_ref
  
  # Build dictionary mapping genes to celltypes
  celltype_marker <- c("Epithelial", "Epithelial", "Macrophages", "Neutrophil", 
                       "mDC", "pDC", "Mast", "T", "NK", "B", "Plasma")
  names(celltype_marker) <- c("TPPP3", "KRT18", "CD68", "FCGR3B", "CD1C",
                              "LILRA4", "TPSB2", "CD3D", "KLRD1", "MS4A1", "IGHG4")
  
  # 
}

###############################################################################

library(rlc)

###############################################################################

# Define the layout, table NxM with N = Rows, M = Cols
openPage( layout = "table8x5", useViewer=FALSE )

###############################################################################

# Define the functions needed in the plots

# 1. Get the average distance between 1000 random marker expressing cells

get_avg_dist = function() {
  a1 <- sample((1:length(subset))[subset_marker], 
               1000, replace=TRUE)
  a2 <- sample((1:length(subset))[subset_marker], 
               1000, replace=TRUE)
  distances <- sqrt(rowSums((pcx[a1,] - pcx[a2])^2))
  return(distances)
}

###############################################################################

# Define the defaults

# Default Marker CD68 for Macrophages
marker <- "CD68"

# Default sample is the first healthy control HC1
sample_select <- "HC3"

# Define the subset
subset <- cells$sample==sample_select

# Default gamma-correction
gamma <- 0.5

# Thus the gamma-correction fraction is
marker_exp <- ((fracs[marker,subset])^gamma)

# Default fraction threshold
threshold <- 0.02

# Thus subset marker is
subset_marker <- (marker_exp > threshold)

# Calculate the nearest neighbors
subset_pcx <<- pcx[subset,]
all_NN <<- RANN::nn2(subset_pcx, k=50)

# Calculate marker NN
marker_NN <<- RANN::nn2(subset_pcx[subset_marker,], k=min(sum(subset_marker), 50))
all_from_marker_NN <<- RANN::nn2(subset_pcx[subset_marker,], query=subset_pcx, k=min(sum(subset_marker), 50))

# Neighbors to consider when looking at the average distance to 10 closest cells
# among the marker expression cells
n_NN <- 25

# Thus the distance between the nearest neighbors is
avg_dist <- rowMeans(marker_NN$nn.dists[,2:(n_NN+1)])
avg_dist_2 <- rowMeans(all_from_marker_NN$nn.dists[,2:(n_NN+1)])

# Default euclidean distance that serves as cutoff to calculate K-Score
cutoff <- 4

# Thus the number of neighbors below the cutoff for each cell is
nn_below <- rowSums(all_NN$nn.dists < cutoff)

# Thus the k-scores are
k_scores <- rowSums(all_from_marker_NN$nn.dists < cutoff)

# Current embedding (default to tsne)

embedding <- tsne

# Define the current k-measure (via global variable)

use_absolut <- TRUE

get_k_measure = function() {
  if (use_absolut) {
    k_measure <<- k_scores
    k_measure_title <<- "Absolut K-Scores"
  } else {
    k_measure <<- k_scores/nn_below
    k_measure_title <<- "Fraction K-Scores"
  }
}

get_k_measure()

# Define above which k-score cells are considered to be part of the class
k_cutoff <- 40

# Thus the classifications are
classification <- (k_measure > k_cutoff)

# Extended classification (the most complicated way possible)
seed_and_classification <- (subset_marker & classification)
only_seed_cells <- (subset_marker & !classification)
only_k_score <- (!subset_marker & classification)

classification_2 <- rep("None", sum(subset))
classification_2[seed_and_classification] <- "Seed+Classification"
classification_2[only_seed_cells] <- "Seed"
classification_2[only_k_score] <- "Classification"

###############################################################################

# Initiate list to stored previous configurations

configurations <- list()

# Create function to update the configuration list

update_configuration  = function() {
  
  # Define name of the list entry
  name <- paste(sample_select, "_", marker, sep="")
  
  # Define empty list which we will add everything too
  entry_list <- list(
    "subset_stored" = subset,
    "gamma_stored" = gamma,
    "marker_exp_stored" = marker_exp,
    "threshold_stored" = threshold,
    "subset_marker_stored" = subset_marker,
    "all_NN_stored" = all_NN,
    "marker_NN_stored" = marker_NN,
    "all_from_marker_NN_stored" = all_from_marker_NN,
    "subset_pcx_stored" = subset_pcx,
    "n_NN_stored" = n_NN,
    "avg_dist_stored" = avg_dist,
    "avg_dist_2_stored" = avg_dist_2,
    "cutoff_stored" = cutoff,
    "nn_below_stored" = nn_below,
    "k_scores_stored" = k_scores,
    "k_cutoff_stored" = k_cutoff,
    "classification_stored" = classification,
    "classification_2_stored" = classification_2
  )
  
  configurations[[name]] <<- entry_list
}

update_configuration()

retrieve_configuraton = function(name) {
  
  # Only called if configuration is available so no sec check included
  
  # Get list entry of interest
  list_oi <- configurations[[name]]
  
  # Define the global variables
  subset <<- list_oi$subset_stored
  gamma <<- list_oi$gamma_stored
  marker_exp <<- list_oi$marker_exp_stored
  threshold <<- list_oi$threshold_stored
  subset_marker <<- list_oi$subset_marker_stored
  all_NN <<- list_oi$all_NN_stored
  marker_NN <<- list_oi$marker_NN_stored
  all_from_marker_NN <<- list_oi$all_from_marker_NN_stored
  subset_pcx <<- list_oi$subset_pcx_stored
  n_NN <<- list_oi$n_NN_stored
  avg_dist <<- list_oi$avg_dist_stored
  avg_dist_2 <<- list_oi$avg_dist_2_stored
  cutoff <<- list_oi$cutoff_stored
  nn_below <<- list_oi$nn_below_stored
  k_scores <<- list_oi$k_scores_stored
  k_cutoff <<- list_oi$k_cutoff_stored
  classification <<- list_oi$classification_stored
  classification_2 <<- list_oi$classification_2_stored
  
  # Update all the charts
  updateCharts()
}

configuration_df = function() {
  configs <- names(configurations)
  vals <- comprehenr::to_df(for (config in configs) configurations[[config]][c(2,4,13,16)])
  samps <- sapply(strsplit(configs, "_"), function(x) x[1])
  genes <- sapply(strsplit(configs, "_"), function(x) x[2])
  df <- as.data.frame(cbind(samps, genes, vals))
  colnames(df) <- c("Sample", "Gene", "Gamma", "Thr", "Cutoff", "K-Cutoff")
  rownames(df) <- c()
  return(df)
}

###############################################################################

# Define update functions

# 1. General update if the marker, the gamma-correction, or the threshold are updated

update_1 = function() {
  
  # Update marker and marker subset
  marker_exp <<- ((fracs[marker,subset])^gamma)
  subset_marker <<- (marker_exp > threshold)
  
  check_empty <- sum(subset_marker)
  
  if (check_empty > n_NN) {
    
    # Calculate marker NN
    marker_NN <<- RANN::nn2(subset_pcx[subset_marker,], k=min(sum(subset_marker), 50))
    all_from_marker_NN <<- RANN::nn2(subset_pcx[subset_marker,], query=subset_pcx, k=min(sum(subset_marker), 50))
    
    # Update k-scores and classification
    avg_dist <<- rowMeans(marker_NN$nn.dists[,2:(n_NN+1)])
    avg_dist_2 <<- rowMeans(all_from_marker_NN$nn.dists[,2:(n_NN+1)])
    nn_below <<- rowSums(all_NN$nn.dists < cutoff)
    k_scores <<- rowSums(all_from_marker_NN$nn.dists < cutoff)
    
    # Get the current k_measure
    get_k_measure()
    
    # Classification
    classification <<- (k_measure > k_cutoff)
    
    # Extended classification (the most complicated way possible)
    seed_and_classification <<- (subset_marker & classification)
    only_seed_cells <<- (subset_marker & !classification)
    only_k_score <<- (!subset_marker & classification)
    
    classification_2 <<- rep("None", sum(subset))
    classification_2[seed_and_classification] <<- "Seed+Classification"
    classification_2[only_seed_cells] <<- "Seed"
    classification_2[only_k_score] <<- "Classification"
    
    #
    display <<- "..."
    
    # Make sure to update the comparison
    get_comparison()
    
    # Update the corresponding charts
    updateCharts()

    } else {
    display <<- "Too few cells expressing the marker (lower the threshold, reduce NN)"
    get_comparison()
    updateCharts(c("display_warning", "hist_frac", "scatter_expression"))
    }
}

# 2. General update if the cutoff or the k_cutoff are changed

update_2 = function() {
  nn_below <<- rowSums(all_NN$nn.dists < cutoff)
  k_scores <<- rowSums(all_from_marker_NN$nn.dists < cutoff)
  
  # Get the current k_measure
  get_k_measure()
  
  # Update classification
  classification <<- (k_measure > k_cutoff)
  
  # Extended classification (the most complicated way possible)
  seed_and_classification <<- (subset_marker & classification)
  only_seed_cells <<- (subset_marker & !classification)
  only_k_score <<- (!subset_marker & classification)
  
  # Update classification
  classification_2 <<- rep("None", sum(subset))
  classification_2[seed_and_classification] <<- "Seed+Classification"
  classification_2[only_seed_cells] <<- "Seed"
  classification_2[only_k_score] <<- "Classification"
  
  # Update the corresponding charts
  updateCharts()
}

###############################################################################

# Define the 4 histograms:

# 1. Histogram of the gamma corrected markergene fractions

lc_hist(
  data = dat(
    value = marker_exp,
    nbins = 100
  ),
  on_clickPosition = function(position) {
    
    # Get only the x position
    x <- position[1]
    
    # Correct for the NA bug
    if (!(is.na(x))) {
      
      # 
      x <- round( x, 3)
     
      # Update the threshold and subset marker
      threshold <<- x
      
      # Update 1
      update_1() 
    }
    
  },
  title = "Gamma-corrected Marker Fraction",
  chartId = "hist_frac",
  transitionDuration = 0,
  place = "A1"
)

lc_vLine(data = dat(
 v = threshold
),
 addLayer = TRUE,
 transitionDuration = 0,
 chartId = "hist_frac"
)

# 2. Histogram to get a feeling for the distances between the marker cells)

lc_hist(
  data = dat(
    value = avg_dist,
    nbins = 100
  ),
  on_clickPosition = function(position) {
    
    # Get only the x position
    x <- position[1]
    
    # Correct for the NA bug
    if (!(is.na(x))) {
      
      # 
      x <- round( x, 1)
      
      # Update the threshold and subset marker
      cutoff <<- x
      
      # Update 2
      update_2()
    }
    
  },
  title = "Average distance from Marker to Marker",
  transitionDuration = 0,
  place = "A2",
  chartId = "hist_M_to_M"
)

lc_vLine(data = dat(
  v = cutoff
),
  addLayer = TRUE,
  chartId = "hist_M_to_M"
)

# 3. Get a feeling for the distances from all cells to the marker cells

lc_hist(
  data = dat(
    value = avg_dist_2,
    nbins = 100
  ),
  on_clickPosition = function(position) {
    
    # Get only the x position
    x <- position[1]
    
    # Correct for the NA bug
    if (!(is.na(x))) {
      
      # 
      x <- round( x, 1)
      
      # Update the threshold and subset marker
      cutoff <<- x
      
      # Update 2
      update_2()
    }
    
  },
  title = "Average distance from All to Marker",
  transitionDuration = 0,
  chartId = "hist_all_to_M",
  place = "B2"
)

lc_vLine(data = dat(
  v = cutoff
),
  addLayer = TRUE,
  chartId = "hist_all_to_M"
)

# 4. Helper Histogram for the number of neighbors which are considered

lc_hist(
  data = dat(
    value = nn_below,
    nbins = 51
  ),
  on_clickPosition = function(position) {
    
    # Get only the x position
    x <- position[1]
    
    # Correct for the NA bug
    if (!(is.na(x))) {
      
      # 
      x <- round( x, 1)
      
      # Update the threshold and subset marker
      k_cutoff <<- x
      
      # Update the k-scores
      update_2()
    }
  },
  title = "Number of Neighbors below Distance Cutoff",
  transitionDuration = 0,
  chartId = "hist_NN_below",
  place = "A3"
  
)

lc_vLine(data = dat(
  v = k_cutoff
),
  addLayer = TRUE,
  transitionDuration = 0,
  chartId = "hist_NN_below"
)

# 5. Histogram for the K-scores

lc_hist(
  data = dat(
    value = k_measure,
    nbins = 51,
    title = k_measure_title
  ),
  on_clickPosition = function(position) {
    
    # Get only the x position
    x <- position[1]
    
    # Correct for the NA bug
    if (!(is.na(x))) {
      
      # 
      x <- round( x, 1)
      
      # Update the threshold and subset marker
      k_cutoff <<- x
      
      # Update 2
      update_2()
    }
  },
  title = "K-Scores",
  transitionDuration = 0,
  chartId = "hist_Kscores",
  place = "A4"
)

lc_vLine(data = dat(
  v = k_cutoff
),
  addLayer = TRUE,
  transitionDuration = 0,
  chartId = "hist_Kscores"
)

###############################################################################

# Add the scatterplot for 2D embeddings

pt_size <- 1

# 1. tSNE with expression of the marker gene

show_comparison <- "expression"

get_comparison <- function() {
  if (show_comparison == "expression") {
    comparison <<- marker_exp
    comparison_title <<- "Gamma-Corrected Marker Fraction"
  }
  else if (show_comparison == "annotation") {
    comparison <<- unlist(cells_ref[subset, "celltype"])
    comparison_title <<-  "Annotation"
  } else if (show_comparison == "library_size") {
    comparison <<- colSums(counts[,subset])
    comparison_title <<-  "Library Size"
  }
}

get_comparison()

lc_scatter(
  dat(
    x = embedding[subset,1],
    y = embedding[subset,2],
    label = colnames(counts)[subset],
    colourValue = comparison,
    title = comparison_title,
    size = pt_size,
    on_marked = function() {
      pass
    }),
  transitionDuration = 0,
  chartId = "scatter_expression",
  place = "A5"
)

# 2. tSNE with classification

lc_scatter(
  dat(
    x = embedding[subset,1],
    y = embedding[subset,2],
    label = colnames(counts)[subset],
    # How to set custom discrete colors in Linked Charts:
    colourValue = classification_2,
    colourDomain = c("Seed+Classification", "Seed", "Classification", "None"),
    palette = c("darkorchid", "blue", "red", "grey"),
    size = pt_size,
    opacity = 0.6,
    on_marked = function() {
      pass
    }),
  transitionDuration = 0,
  title = "Classification",
  chartId = "scatter_classification",
  place = "B5"
)

# 3. tSNE with number of neighbors considering the current cutoff
lc_scatter(
  dat(
    x = embedding[subset,1],
    y = embedding[subset,2],
    label = colnames(counts)[subset],
    colourValue = nn_below,
    size = pt_size,
    on_marked = function() {
      pass
    }),
  transitionDuration = 0,
  title = "Number of Neigbors below Distance Cutoff",
  chartId = "scatter_NN_below",
  place = "B3"
)

# 4. tSNE with the k-score as color
lc_scatter(
  dat(
    x = embedding[subset,1],
    y = embedding[subset,2],
    label = colnames(counts)[subset],
    colourValue = k_measure,
    size = pt_size,
    title = k_measure_title,
    on_marked = function() {
      pass
    }),
  transitionDuration = 0,
  title = "K-scores",
  chartId = "scatter_Kscores",
  place = "B4"
)


###############################################################################

# Define the input functionality:

# Define the width of the input elements
elem_width <- 300

# 1. Which sample is used (every chart will have to be updated)

lc_input(type = "text", label=c(""),
         dat(value = sample_select),
         on_click = function(smp) {
           
           # Add security check, whether the gene really exists
           if (!(smp %in% unique(cells$sample))) {
             display <<- "Sample does not exist"
             updateCharts()
             return()
           }
           
           # Store the configuration
           update_configuration()
           
           # Update selected sample
           sample_select <<- smp
           
           # Check if configuration is already available
           name <- paste(sample_select, "_", marker, sep="")
           
           # Get the configuration names
           configs_available <- names(configurations)
           splitted <- strsplit(configs_available, "_")
           smp_check <- grep(sample_select, splitted)[1]
          
           # If possible get the whole configuration
           if (name %in% names(configurations)) {
             
             # Retrieve Configuration
             retrieve_configuraton(name)
             
           # If the sample has been selected beofre get, the sample specific things
           } else if (!(is.na(smp_check))) {
             
             # Retrieve subsets and all_NN
             subset <<- configurations[[smp_check]]$subset_stored
             subset_pcx <<- configurations[[smp_check]]$subset_pcx_stored
             all_NN <<- configurations[[smp_check]]$all_NN_stored
               
               
             # Calculate everything else -> Update 1
             update_1()
           
           # Otherwise calculate everything
           } else {
             
             # Update subset
             subset <<- cells$sample==sample_select
             
             # Calculate the nearest neighbors
             subset_pcx <<- pcx[subset,]
             all_NN <<- RANN::nn2(subset_pcx, k=50)
             
             # Update 1
             update_1()
           }
         },
         place = "C1",
         chartId = "input_sample",
         width = elem_width,
         title = "Sample"
)

# 2. Which marker gene is used (every chart will have to be updated)

lc_input(type = "text", label=c(""),
         dat(value = marker),
         on_click = function(gene) {
           
           # Add security check, whether the gene really exists
           if (!(gene %in% rownames(fracs))) {
             display <<- "Gene does not exist"
             updateCharts(c("display_warning"))
             return()
           }
           
           # Store the configuration
           update_configuration()
           
           # Update marker
           marker <<- gene
           
           # Check if configuration is already available
           name <- paste(sample_select, "_", marker, sep="")
           if (name %in% names(configurations)) {
             # Retrieve Configuration
             retrieve_configuraton(name)
            } else {
           
             # Update 1
             update_1()
            }
          }, 
         place = "D1",
         chartId = "input_marker",
         width = elem_width,
         title = "Marker Gene"
)

# 3. Which gamma correction factor should be used (every chart will have to be updated)

lc_input(type = "range",
         on_click = function(g) {
           
           # Update gamma
           gamma <<- g
           
           # Update 1
           update_1()
           
         },
         dat(value = gamma), 
         min = c(0), 
         max = c(1),
         step=c(0.1),
         chartId = "input_gamma",
         place = "C2",
         width = elem_width,
         title = "Gamma",
         labels = c("")
)

# 4. Which threshold to use to consider cells as marker expressing ones

lc_input(type = "text", 
         data = dat(value = threshold),
         on_click = function(thr) {
           
            # Update the threshold
            threshold <<- as.numeric(thr)
           
             # Update 1
            update_1()
           
            },
         chartId = "input_threshold",
         place = "D2",
         width = elem_width,
         title = "Threshold",
         labels = c("")
)

# 5. Which distance to use to show the average distances

lc_input(type = "text", 
         dat(value = n_NN),
         on_click = function(n_neighbors) {
           
           # Update how many neighbors to consider
           n_NN <<- as.numeric(n_neighbors)
           
           # Update k-scores and classifcation
           avg_dist <<- rowMeans(marker_NN$nn.dists[,2:(n_NN+1)])
           avg_dist_2 <<- rowMeans(all_from_marker_NN$nn.dists[,2:(n_NN+1)])
           
           # Update display to be sure
           display <<- "..."
           
           # Update the corresponding charts
           updateCharts()
           
         }, 
         place = "C3",
         chartId = "input_nn_cons",
         width = elem_width,
         title = "Considered Neighbors",
         labels = c("")
)

# 6. Which distance cutoff to use to calculate the k-score

lc_input(type = "text", 
         dat(value = cutoff),
         on_click = function(cutoff) {
           
           # Update the cutoff
           cutoff <<- as.numeric(cutoff)
           
          # Update 2
           update_2()

         }, 
         place = "D3",
         chartId = "input_cutoff",
         width = elem_width,
         title = "Distance Cutoff",
         labels = c("")
)

# 7. Which k-cutoff to use for classification

lc_input(type = "text", 
         dat(value = k_cutoff),
         on_click = function(entered) {
           k_cutoff <<- as.numeric(entered)
            
           # Update the k-scores
           classification <<- (k_scores > k_cutoff)
           
           # Extended classification (the most complicated way possible)
           seed_and_classification <<- (subset_marker & classification)
           only_seed_cells <<- (subset_marker & !classification)
           only_k_score <<- (!subset_marker & classification)
           
           classification_2 <<- rep("None", sum(subset))
           
           for (i in 1:sum(subset)) {
             if (seed_and_classification[i]) {
               classification_2[i] <<- "Seed+Classification"
             } else if (only_seed_cells[i]) {
               classification_2[i] <<- "Seed"
             } else if (only_k_score[i]) {
               classification_2[i] <<- "Classification"
             }
           }
           
           # Update the corresponding charts
           updateCharts()
           
         }, 
         place = "C4",
         chartId = "input_Kcutoff",
         width = elem_width,
         title = "K-Score Cutoff",
         labels = c("")
)

# 8. Size of points in the tSNE plots

lc_input(type = "range",
         on_click = function(pt) {
           
           # Update point size
           pt_size <<- pt
           
           # Update the tSNEs
           updateCharts(c("scatter_expression", "scatter_classification", 
                          "scatter_NN_below", "scatter_Kscores"))
           
         },
         dat(value = pt_size), 
         min = c(1), 
         max = c(6),
         step=c(1),
         chartId = "input_pt",
         place = "D4",
         width = elem_width,
         title = "Point Size",
         labels = c("")
)

# 9. Change between UMAP and tSNE

lc_input(type = "radio", labels = c("t-SNE", "UMAP"), 
         on_click = function(value) {
           
           # Check which button was clicked and change to the corresponding embedding
           if (value == 1) {
             
             embedding <<- tsne
             
           } else if (value == 2) {
             
             embedding <<- umap
             
           }
           
           # Update the scatterplots
           updateCharts(c("scatter_expression", "scatter_classification", 
                        "scatter_NN_below", "scatter_Kscores"))
           
         },
         value = 1, 
         place = "C5",
         title = "2D-Embedding"
         )

# 10. Button to load sleepwalker

lc_input(type = "button", labels = c("Load Sleepwalker"), 
         on_click = function(value) {
           sleepwalk::sleepwalk(embedding[subset,], pcx[subset,])
         },
         place = "D5")

# 11. Change between k-scores and k-score fraction

lc_input(type = "radio", labels = c("Absolut K-Scores", "Fraction K-Scores"), 
         on_click = function(value) {
           
           # Check which button was clicked and change to the corresponding embedding
           if (value == 1) {
             
             use_absolut <<- TRUE
             
           } else if (value == 2) {
             
             use_absolut <<- FALSE
           }
           
           # Get the updated k_measure
           get_k_measure()
           
           # Update the classification
           classification <<- (k_measure > k_cutoff)
           
           # Extended classification (the most complicated way possible)
           seed_and_classification <<- (subset_marker & classification)
           only_seed_cells <<- (subset_marker & !classification)
           only_k_score <<- (!subset_marker & classification)
           
           classification_2 <<- rep("None", sum(subset))
           classification_2[seed_and_classification] <<- "Seed+Classification"
           classification_2[only_seed_cells] <<- "Seed"
           classification_2[only_k_score] <<- "Classification"
           
           # Update the scatterplots
           updateCharts()
           
         },
         value = 1, 
         place = "E4",
         title = "K-Measure")

# 12. Change between Expression, Annotation, and Library Size

lc_input(type = "radio", labels = c("Expression", "Annotation", "Library Size"), 
         on_click = function(value) {
           
           # Check which button was clicked and change to the corresponding embedding
           if (value == 1) {

             show_comparison <<- "expression"
             
           } else if (value == 2) {
             
             show_comparison <<- "annotation" 
             
           } else if (value == 3) {
             
             show_comparison <<- "library_size"
               
           }
           
           # Update the comparison
           get_comparison()
           
           # Update the scatterplots
           updateCharts(c("scatter_expression"))
           
         },
         value = 1, 
         place = "E5",
         title = "Comparison",
         labels = c("")
         )

###############################################################################


# Include display for test purposes

display <- "..."

lc_html(dat(
  content = data.frame(Warning=display)), 
  chartId = "display_warning",
  place = "E1"
)

# Include second display to show the classification data

lc_html(dat(
  content = data.frame("Seed_Class" = sum(seed_and_classification),
                       "Only_Seed" = sum(only_seed_cells),
                       "Sum_Seed" = sum(subset_marker),
                       "Only_Class" = sum(only_k_score),
                       "Sum_Class" = sum(classification),
                       "None" = sum(classification_2 == "None"),
                       "Total" = length(classification_2))), 
  place =  "E2"
)

# Include third display to show stored configurations

lc_html(dat(
  content = configuration_df()
  ), 
  place =  "F1"
)
