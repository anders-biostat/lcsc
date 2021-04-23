
###############################################################################

# Load library

library(sparseMatrixStats)
library(Matrix)

###############################################################################

# Load Data

# Specify input directory
dir <- "~/3_r_projects/Anders_Lab_Local/22-04-21_processing/data/"  # windows

load_flag = FALSE

if (load_flag) {
  load(paste(dir,"counts.RData", sep=""))
  load(paste(dir,"fracs.RData", sep=""))
  load(paste(dir,"cells.RData", sep=""))
  load(paste(dir,"sample_info.RData", sep=""))
  load(paste(dir,"pcx.RData", sep=""))
  load(paste(dir,"pcs.RData", sep=""))
  load(paste(dir,"tsne_pcx.RData", sep=""))
  load(paste(dir,"tsne_pcs.RData", sep=""))
  load(paste(dir,"umap_pcx.RData", sep=""))
  load(paste(dir,"umap_pcs.RData", sep=""))
  load(paste(dir,"all_NN_pcx.RData", sep=""))
  load(paste(dir,"all_NN_pcs.RData", sep=""))
  load(paste(dir,"cells_ref.RData", sep=""))
  load(paste(dir,"celltype_marker.RData", sep=""))
}

###############################################################################

library(rlc)

###############################################################################

# Define the layout, table NxM with N = Rows, M = Cols
openPage( layout = "table8x6", useViewer=FALSE )

###############################################################################

# Define the defaults

# Default Marker CD68 for Macrophages
marker <- "CD68"

# Default sample is the first healthy control HC1
sample_select <- "HC3"

# Define the subset
subset <- cells$sample==sample_select

# Which PC space to use
pc_space <- pcx
all_NN <- all_NN_pcx

# Getting the nearest neighbors
correction <- sample_info$cumulative[grep(sample_select, sample_info$sample)]
all_NN_subset <- list(idx = (all_NN$idx[subset,] - correction), 
                      dist = all_NN$dist[subset,])

# Get the average distance for cell_i to x nearest neighbors
n_NN <- 25
avg_dist <- rowMeans(all_NN_subset$dist[,2:(n_NN+1)])

# Define the distance cutoff
cutoff <- 6

# Calculate number of neighbors below cutoff
nn_below <- rowSums(all_NN_subset$dist < cutoff)

# Get the counts for the marker genes for all cells from the sample
c_M <- counts[marker, subset]

# For each cell of the sample, get the sum of the counts for the marker 
# genes for cells that are closer than the cutoff
counts_ci <- sapply(1:nrow(all_NN_subset$idx), function(cell_i){
  # Get the indices of all neighbors that are closer than the cutoff
  # (including cell_i !)
  below_cutoff <- all_NN_subset$idx[cell_i, all_NN_subset$dist[cell_i,] < cutoff]
  # Get the sum of counts
  sum(c_M[below_cutoff])
})

# Get the total counts for all the cells in the sample
c_Tot <- colSums(counts[, subset])

# For each cell of the sample, get the total counts
# for cells that are closer than the cutoff
totals_ci <- sapply(1:nrow(all_NN_subset$idx), function(cell_i){
  # Get the indices of all neighbors that are closer than the cutoff
  # (including cell_i !)
  below_cutoff <- all_NN_subset$idx[cell_i, all_NN_subset$dist[cell_i,] < cutoff]
  # Get the sum of counts
  sum(c_Tot[below_cutoff])
})

# Calculate the fractions
f_M <- counts_ci/totals_ci

# Default gamma-correction
gamma <- 0.5

# Calculate the corrected ratio
f_M_c <- f_M^gamma

# Default fraction threshold
threshold <- 0.02

# Classification
classification <- (f_M_c > threshold)

###############################################################################

# Initial constant for figures

# Current embedding (default to tsne)
embedding_str <- "tsne"
tsne <- tsne_pcx
embedding <- tsne

# Initial point size
pt_size <- 1

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
    "threshold_stored" = threshold,
    "all_NN_subset_stored" = all_NN_subset,
    "n_NN_stored" = n_NN,
    "avg_dist_stored" = avg_dist,
    "counts_ci_stored" = counts_ci,
    "totals_ci_stored" = totals_ci,
    "f_M_stored" = f_M,
    "f_M_c_stored" = f_M_c,
    "cutoff_stored" = cutoff,
    "classification_stored" = classification
  )
  
  configurations[[name]] <<- entry_list
}

update_configuration()

retrieve_configuraton = function(name) {
  
  # Get list entry of interest
  list_oi <- configurations[[name]]
  
  # Define the global variables
  subset <<- list_oi$subset_stored
  gamma <<- list_oi$gamma_stored
  threshold <<- list_oi$threshold_stored
  all_NN_subset <<- list_oi$all_NN_subset_stored
  n_NN <<- list_oi$n_NN_stored
  avg_dist <<- list_oi$avg_dist_stored
  counts_ci <<- list_oi$counts_ci_stored
  totals_ci_stored <<- list_oi$totals_ci_stored
  f_M <<- list_oi$f_M_stored
  f_M_c <<- list_oi$f_M_c_stored
  cutoff <<- list_oi$cutoff_stored
  classification <<- list_oi$classification_stored
  
  # Update all the charts
  updateCharts()
}

# Define function to store the settings for each sample-gene combination
configuration_df = function() {
  configs <- names(configurations)
  vals <- comprehenr::to_df(for (config in configs) configurations[[config]][c(2,3,11)])
  samps <- sapply(strsplit(configs, "_"), function(x) x[1])
  genes <- sapply(strsplit(configs, "_"), function(x) x[2])
  df <- as.data.frame(cbind(samps, genes, vals))
  colnames(df) <- c("Sample", "Gene", "Gamma", "Threshold", "Cutoff")
  rownames(df) <- c()
  return(df)
}

###############################################################################

# 1. Histogram to get feeling for the distances between the cells

lc_hist(
  data = dat(
    value = avg_dist,
    nbins = 100,
    title = paste("Average distance to the", n_NN, "nearest neighbors")
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
      
      # Update 1
      update_1()
    }
    
  },
  transitionDuration = 0,
  place = "A1",
  chartId = "hist_dists"
)

lc_vLine(data = dat(
  v = cutoff
),
addLayer = TRUE,
chartId = "hist_dists"
)

# 2. Histogram to show how many neighbors each cell has below the cutoff

lc_hist(
  data = dat(
    value = nn_below,
    nbins = 100
  ),
  title = "Number of neighbors below the distance cutoff",
  transitionDuration = 0,
  place = "A2",
  chartId = "hist_nn_below"
)

# 3. 2D DimPlot to how many neighbors each cell has below the cutoff

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
  chartId = "dimred_NN_below",
  place = "B2"
)

# 4. Histogram to show the fraction

lc_hist(
  data = dat(
    value = f_M_c,
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
      
      # Update the classification
      classification <<- (f_M_c > threshold)
      
      # Update charts
      updateCharts()
    }
    
  },
  title = "Gamma-corrected Fraction",
  chartId = "hist_frac",
  transitionDuration = 0,
  place = "A3"
)

lc_vLine(data = dat(
  v = threshold
),
addLayer = TRUE,
transitionDuration = 0,
chartId = "hist_frac"
)

# 5. Scatterplot to show the fractions

lc_scatter(
  data = dat(
    x = totals_ci^gamma,
    y = f_M_c,
    size = pt_size
  ),
  on_clickPosition = function(position) {
    
    # Correct for the NA bug
    if (!(is.na(position[1]))) {
      
      # 
      x <- round( position[2], 3)
      
      # Update the threshold and subset marker
      threshold <<- x
      
      # Update the classification
      classification <<- (f_M_c > threshold)
      
      # Update charts
      updateCharts()
    }
  },
  transitionDuration = 0,
  title = "Fraction vs. Library Size",
  chartId = "scatter_frac",
  place = "B3"
  # TODO: Add x and y labels, check whether is right
)

lc_hLine(data = dat(
  h = threshold
  ),
  addLayer = TRUE,
  transitionDuration = 0,
  chartId = "scatter_frac"
)

# 6. 2D DimPlot to how show the fraction

lc_scatter(
  dat(
    x = embedding[subset,1],
    y = embedding[subset,2],
    label = colnames(counts)[subset],
    colourValue = f_M_c,
    size = pt_size,
    on_marked = function() {
      pass
    }),
  transitionDuration = 0,
  title = "Fraction",
  chartId = "dimred_frac",
  place = "A4"
)

# 7. 2D DimPlot to how show the classification

lc_scatter(
  dat(
    x = embedding[subset,1],
    y = embedding[subset,2],
    label = colnames(counts)[subset],
    colourValue = classification,
    size = pt_size,
    on_marked = function() {
      pass
    }),
  transitionDuration = 0,
  title = "Classification",
  chartId = "dimred_class",
  place = "B4"
)

# 8. 2D DimPlot to show the counts

lc_scatter(
  dat(
    x = embedding[subset,1],
    y = embedding[subset,2],
    label = colnames(counts)[subset],
    colourValue = counts[marker,subset],
    size = pt_size,
    on_marked = function() {
      pass
    }),
  transitionDuration = 0,
  title = "Counts",
  chartId = "dimred_counts",
  place = "A5"
)

# 9. 2D DimPlot to show library size

lc_scatter(
  dat(
    x = embedding[subset,1],
    y = embedding[subset,2],
    label = colnames(counts)[subset],
    colourValue = colSums(counts[,subset]),
    size = pt_size,
    on_marked = function() {
      pass
    }),
  transitionDuration = 0,
  title = "Library Size (total UMI count)",
  chartId = "dimred_size",
  place = "B5"
)

# 10. 2D DimPlot to show the annotation

lc_scatter(
  dat(
    x = embedding[subset,1],
    y = embedding[subset,2],
    label = colnames(counts)[subset],
    colourValue = cells_ref$celltype[subset],
    size = pt_size,
    on_marked = function() {
      pass
    }),
  transitionDuration = 0,
  title = "Annotation",
  chartId = "dimred_annotation",
  place = "A6"
)

###############################################################################

# Define the input

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
           
           # Store the previous configuration
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
             
             # Otherwise make everything new
           } else {
             
             # Define the subset
             subset <<- cells$sample==sample_select
             
             # Getting the nearest neighbors
             correction <<- sample_info$cumulative[grep(sample_select, sample_info$sample)]
             all_NN_subset <<- list(idx = (all_NN$idx[subset,] - correction), 
                                   dist = all_NN$dist[subset,])
             
             # Get the average distance for cell_i to x nearest neighbors
             avg_dist <<- rowMeans(all_NN_subset$dist[,2:(n_NN+1)])
             
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

# 3. Which distance to use to show the average distances

lc_input(type = "text", 
         dat(value = n_NN),
         on_click = function(n_neighbors) {
           
           if (!(is.numeric(n_neighbors))) {
             display <<- "Enter a number"
             updateCharts(c("display_warning"))
             return()
           }
           
           # Checking
           if (n_neighbors >= 100) {
             display <<- "Neighbors considered must be below 100"
             updateCharts(c("display_warning"))
             return()
           }
           
           # Update how many neighbors to consider
           n_NN <<- as.numeric(n_neighbors)
           
           # Update k-scores and classification
           avg_dist <<- rowMeans(all_NN_subset$dist[,2:(n_NN+1)])
           
           # Update display to be sure
           display <<- "..."
           
           # Update the corresponding charts
           updateCharts()
           
         }, 
         place = "C2",
         chartId = "input_nn_cons",
         width = elem_width,
         title = "Considered Neighbors",
         labels = c("")
)

# 4. Which distance cutoff to use to calculate the k-score

lc_input(type = "text", 
         dat(value = cutoff),
         on_click = function(cutoff) {
           
           if (!(is.numeric(cutoff))) {
             display <<- "Enter a number"
             updateCharts(c("display_warning"))
             return()
           }
           
           # Update the cutoff
           cutoff <<- as.numeric(cutoff)
           
           # Update 2
           update_1()
           
         }, 
         place = "D2",
         chartId = "input_cutoff",
         width = elem_width,
         title = "Distance Cutoff",
         labels = c("")
)

# 5. Which gamma correction factor should be used (every chart will have to be updated)

lc_input(type = "range",
         on_click = function(g) {
           
           # Update gamma
           gamma <<- g
           
           # Calculate the corrected fractions
           f_M_c <<- f_M^gamma
           
           # Update the classification
           classification <<- (f_M_c > threshold)
           
           # Update charts
           updateCharts()
           
         },
         dat(value = gamma), 
         min = c(0), 
         max = c(1),
         step=c(0.1),
         chartId = "input_gamma",
         place = "C3",
         width = elem_width,
         title = "Gamma",
         labels = c("")
)

# 6. Which threshold to use to consider cells as marker expressing ones

lc_input(type = "text", 
         data = dat(value = threshold),
         on_click = function(thr) {
           
           if (!(is.numeric(thr))) {
             display <<- "Enter a number"
             updateCharts(c("display_warning"))
             return()
           }
           
           # Update the threshold
           threshold <<- thr
           
           # Update the classification
           classification <<- (f_M_c > threshold)
           
           # Update charts
           updateCharts()
           
         },
         chartId = "input_threshold",
         place = "D3",
         width = elem_width,
         title = "Threshold",
         labels = c("")
)

# 7. Size of points in the scatter plots

lc_input(type = "range",
         on_click = function(pt) {
           
           # Update point size
           pt_size <<- pt
           
           # Update the pt size
           updateCharts()
           
         },
         dat(value = pt_size), 
         min = c(1), 
         max = c(6),
         step=c(1),
         chartId = "input_pt",
         place = "C4",
         width = elem_width,
         title = "Point Size",
         labels = c("")
)

# 8. Change between UMAP and tSNE

lc_input(type = "radio", labels = c("t-SNE", "UMAP"), 
         on_click = function(value) {
           
           # Check which button was clicked and change to the corresponding embedding
           if (value == 1) {
             
             embedding_str <<- "tsne" # little helper
             embedding <<- tsne
             
           } else if (value == 2) {
             
             embedding_str <<- "umap" # little helper
             embedding <<- umap
             
           }
           
           # Update the scatterplots
           updateCharts()
           
         },
         value = 1, 
         place = "D4",
         title = "2D-Embedding"
)

# 9. Change between PCX and PCS space
# TODO
lc_input(type = "radio", labels = c("PCX (mutual)", "PCS (sample specific)"), 
         on_click = function(value) {
           
           # Check which button was clicked and change to the corresponding embedding
           if (value == 1) {
             
             # Which PC space to use
             pc_space <<- pcx
             all_NN <<- all_NN_pcx
             tsne <<- tsne_pcx
             umap <<- umap_pcx
             
             # Update the embedding as well
             ifelse(embedding_str == "tsne", embedding <<- tsne, embedding <<- umap)
             
             # Getting the nearest neighbors
             correction <<- sample_info$cumulative[grep(sample_select, sample_info$sample)]
             all_NN_subset <<- list(idx = (all_NN$idx[subset,] - correction), 
                                   dist = all_NN$dist[subset,])
             
             # Recalculate avg_dist
             avg_dist <<- rowMeans(all_NN_subset$dist[,2:(n_NN+1)])
             
             # Update everything else
             update_1()
             
           } else if (value == 2) {
             
             # Which PC space to use
             pc_space <<- pcs
             all_NN <<- all_NN_pcs
             tsne <<- tsne_pcs
             umap <<- umap_pcs
             
             # Update the embedding as well
             ifelse(embedding_str == "tsne", embedding <<- tsne, embedding <<- umap)
             
             # Getting the nearest neighbors
             correction <<- sample_info$cumulative[grep(sample_select, sample_info$sample)]
             all_NN_subset <<- list(idx = (all_NN$idx[subset,] - correction), 
                                    dist = all_NN$dist[subset,])
             
             # Recalculate avg_dist
             avg_dist <<- rowMeans(all_NN_subset$dist[,2:(n_NN+1)])
             
             # Update everything else
             update_1()
             
           }
           
           # Update the scatterplots
           updateCharts()
           
         },
         value = 1, 
         place = "D5",
         title = "PC Space"
)

# 10. Button to load sleepwalker

lc_input(type = "button", labels = c("Load Sleepwalker"), 
         on_click = function(value) {
           sleepwalk::sleepwalk(embedding[subset,], pc_space[subset,])
         },
         place = "C5")

###############################################################################

# Define update functions

# Update 1, which includes everything post sample selection
# Will need to be called if the whole configuration is not already available

update_1 = function() {
  
  # Calculate number of neighbors below cutoff
  nn_below <<- rowSums(all_NN_subset$dist < cutoff)
  
  # Get the counts for the marker genes for all cells from the sample
  c_M <<- counts[marker, subset]
  
  # For each cell of the sample, get the sum of the counts for the marker 
  # genes for cells that are closer than the cutoff
  counts_ci <<- sapply(1:nrow(all_NN_subset$idx), function(cell_i){
    # Get the indices of all neighbors that are closer than the cutoff
    # (including cell_i !)
    below_cutoff <- all_NN_subset$idx[cell_i, all_NN_subset$dist[cell_i,] < cutoff]
    # Get the sum of counts
    sum(c_M[below_cutoff])
  })
  
  # Get the total counts for all the cells in the sample
  c_Tot <<- colSums(counts[, subset])
  
  # For each cell of the sample, get the total counts
  # for cells that are closer than the cutoff
  totals_ci <<- sapply(1:nrow(all_NN_subset$idx), function(cell_i){
    # Get the indices of all neighbors that are closer than the cutoff
    # (including cell_i !)
    below_cutoff <- all_NN_subset$idx[cell_i, all_NN_subset$dist[cell_i,] < cutoff]
    # Get the sum of counts
    sum(c_Tot[below_cutoff])
  })
  
  # Calculate the fractions
  f_M <<- counts_ci/totals_ci
  
  # Calculate the corrected fractions
  f_M_c <<- f_M^gamma
  
  # Update the classification
  classification <<- (f_M_c > threshold)
  
  # Update the display
  display <<- "..."
  
  # Update charts
  updateCharts()
}

###############################################################################

# Add html output

display <- "..."

lc_html(dat(
  content = data.frame(Warning=display)), 
  chartId = "display_warning",
  place = "E1"
)


# Include display to show stored configurations

lc_html(dat(
  content = configuration_df()
), 
place =  "F1"
)

# Include display to show classification

lc_html(dat(
  content = data.frame(Positive=(sum(classification)),
                       Negative=sum(!(classification)))
), 
place =  "E2"
)









