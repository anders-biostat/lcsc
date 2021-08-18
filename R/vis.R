
#' @importFrom dplyr %>%
#' @export
dplyr::`%>%`

#' Visualize and annotate single cell data
#' @name lc_vis
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @param cells data.frame with cell metadata
#' @param counts count matrix (cols = cells, genes = rows)
#' @param pc_space PCA coordinates
#' @param embedding coordinates for two-dimensional embedding
#' @param nn list containing nearest neighbor information
#' @param k number of neighbors considered for smoothing
#' @export
lc_vis <- function(cells, counts, pc_space, embedding, nn, k=50){

  ###############################################################################

  # Generate sample info
  tibble::tibble(sample = unique(cells$sample)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(total = sum(cells$sample == sample)) %>%
    dplyr::ungroup() -> sample_info

  # Correction factor adjust indices for nn below
  sample_info$correction = c(0, cumsum(sample_info$total)[1:nrow(sample_info)-1])

  ###############################################################################

  # Start rlc page
  startPage = system.file("http_root/rlc_layout.html", package="lcsc")
  rlc::openPage(startPage=startPage, useViewer=FALSE)

  ###############################################################################

  # Initialize
  g <- globalenv()
  marker <- sample(rownames(counts),1)    # select random gene
  sample_select <- sample_info$sample[1]
  selection <- cells$sample==sample_select

  # Getting the nearest neighbors of all cells in the subset.
  correction <- sample_info %>% dplyr::filter(sample==sample_select) %>% dplyr::pull(correction)
  nn_subset <- list(idx = (nn$idx[selection,] - correction),
                    dists = nn$dist[selection,])

  # Get the average distance for cell_i to x nearest neighbors.
  n_NN <- 25
  avg_dist <- Matrix::rowMeans(nn_subset$dists[,2:(n_NN+1)])

  # Calculate marker counts for cell i also considering the neighboring cells
  c_M <- counts[marker, selection]
  counts_ci <- sapply(1:nrow(nn_subset$idx), function(cell_i){
    sum(c_M[nn_subset$idx[cell_i,]])
  })

  # Calculate the total counts for cell i also condering the neigboring cells
  c_tot <- Matrix::colSums(counts[, selection])
  totals_ci <- sapply(1:nrow(nn_subset$idx), function(cell_i){
    sum(c_tot[nn_subset$idx[cell_i,]])
  })

  # Calculate the gamma corrected fractions
  f_M <- counts_ci/totals_ci
  gamma = 0.5
  f_M_c <- f_M^gamma

  # Default fraction threshold
  threshold <- 0.02

  # Default positive classification if above threshold
  above <- TRUE

  # Classification
  ifelse(above, classification <- f_M_c > threshold, classification <- f_M_c < threshold)

  ###############################################################################

  # Initial constant for figures

  # Point size
  pt_size <- 1

  ###############################################################################

  # Initialize the plots

  # 1. Histogram to get feeling for the distances between the cells

  rlc::lc_hist(
    titleSize=16,
    width=320,
    height=320,
    data = rlc::dat(
      value = avg_dist,
      nbins = 100,
      title = paste("Mean dist to", n_NN, "th nearest neigbor")
    ),
    transitionDuration = 0,
    place = "avg_dist_hist",
    chartId = "hist_dists"
  )

  # 2. Embedding to show the average distance of each cell

  rlc::lc_scatter(
    titleSize=16,
    width=320,
    height=320,
    rlc::dat(
      x = embedding[selection,1],
      y = embedding[selection,2],
      label = colnames(counts)[selection],
      colourValue = avg_dist,
      size = pt_size,
      title = paste("Mean dist to", n_NN, "th nearest neigbor"),
      on_marked = function() {
      }),
    transitionDuration = 0,
    chartId = "nn_below_plot",
    place = "nn_below_plot"
  )

  # 3. Plots to show the distribution of the marker gene fraction
  # Default fraction is scatter

  rlc::lc_scatter(
    titleSize=16,
    width=320,
    height=320,
    data = rlc::dat(
      x = totals_ci^gamma,
      y = f_M_c,
      size = pt_size,
      colourValue = classification,
      colourDomain = c(TRUE, FALSE),
      palette = c("green", "red")
    ),
    on_clickPosition = function(position) {

      # Correct for the NA bug (on_ClickPosition return first NA and then the value)
      if (!(is.na(position[1]))) {

        # Round to reduce number of digit after decimal point
        x <- position[2]

        # Update the threshold
        threshold <<- x

        # Update the classification
        ifelse(above, classification <<- f_M_c > threshold, classification <<- f_M_c < threshold)

        # Update charts
        rlc::updateCharts()
      }
    },
    transitionDuration = 0,
    title = "Fraction vs. Library Size",
    chartId = "fraction_plot",
    place = "fraction_plot",
    on_marked = function() {
      rlc::mark(rlc::getMarked(chartId="fraction_plot", layerId = "Layer1"), "classification_plot", clear=TRUE)
    }
  )

  rlc::lc_hLine(data = rlc::dat(
    h = threshold
  ),
  addLayer = TRUE,
  transitionDuration = 0,
  chartId = "fraction_plot"
  )

  # 4. 2D DimPlot to how show the fraction

  rlc::lc_scatter(
    titleSize=16,
    width=320,
    height=320,
    rlc::dat(
      x = embedding[selection,1],
      y = embedding[selection,2],
      label = colnames(counts)[selection],
      colourValue = f_M_c,
      size = pt_size,
      on_marked = function() {
      }),
    transitionDuration = 0,
    title = "Fraction",
    chartId = "dimred_frac",
    place = "fraction_embedding"
  )

  # 5. Plots to show the classification
  # (default plot is showing classification based on current rule)

  rlc::lc_scatter(
    titleSize=16,
    width=320,
    height=320,
    rlc::dat(
      x = embedding[selection,1],
      y = embedding[selection,2],
      label = colnames(counts)[selection],
      colourValue = classification,
      colourDomain = c(TRUE, FALSE),
      palette = c("green", "red"),
      size = pt_size,
      on_marked = function() {
      }),
    transitionDuration = 0,
    title = "Classification Rule",
    chartId = "classification_plot",
    place = "classification_plot"
  )

  # 6. Plots to show additional information in the 2D-embedding
  # (default is the number of count for the marker gene != smoothed fraction)

  rlc::lc_scatter(
    titleSize=16,
    width=320,
    height=320,
    rlc::dat(
      x = embedding[selection,1],
      y = embedding[selection,2],
      label = colnames(counts)[selection],
      colourValue = counts[marker,selection],
      size = pt_size,
      on_marked = function() {
      }),
    transitionDuration = 0,
    title = "Counts",
    chartId = "additional_data_plot",
    place = "additional_data_plot"
  )

  ###############################################################################

  # Define the input

  # 1. Which sample is used.

  rlc::lc_input(type = "text", label=c(""),
                rlc::dat(value = sample_select),
                on_click = function(smp) {

                  # Add security check, whether the gene really exists
                  if (!(smp %in% sample_info$sample)) {
                    display <<- "Sample does not exist"
                    rlc::updateCharts()
                    return()
                  }

                  # Update selected sample
                  sample_select <<- smp
                  selection <<- cells$sample==sample_select

                  # Getting the nearest neighbors of all cells in the subset.
                  correction <<- sample_info %>% dplyr::filter(sample==sample_select) %>% dplyr::pull(correction)
                  nn_subset <<- list(idx = (nn$idx[selection,] - correction),
                                     dists = nn$dist[selection,])

                  # Get the average distance for cell_i to x nearest neighbors
                  avg_dist <<- Matrix::rowMeans(nn_subset$dist[,2:(n_NN+1)])

                  # Update
                  update_1()


                },
                place = "sample",
                chartId = "input_sample",
                width = 240,
                title = "Sample"
  )

  # 2. Which marker gene is used.

  rlc::lc_input(type = "text", label=c(""),
                rlc::dat(value = marker),
                on_click = function(gene) {

                  # Add security check, whether the gene really exists
                  if (!(gene %in% rownames(counts))) {
                    display <<- "Gene does not exist"
                    rlc::updateCharts(c("display_warning"))
                    return()
                  }

                  # Update marker
                  marker <<- gene

                  # Update
                  update_1()

                },
                place = "gene",
                chartId = "input_marker",
                width = 240,
                title = "Marker Gene"
  )

  # 3. How many neighbors to consider to show average distance

  rlc::lc_input(type = "text",
                rlc::dat(value = n_NN),
                on_click = function(n_neighbors) {

                  if (!(is.numeric(n_neighbors))) {
                    display <<- "Enter a number"
                    rlc::updateCharts(c("display_warning"))
                    return()
                  }

                  # Checking I
                  if (n_neighbors >= k) {
                    display <<- paste("Neighbors considered must be below", k)
                    rlc::updateCharts(c("display_warning"))
                    return()
                  }

                  # Checking II
                  if (n_neighbors < 2) {
                    display <<- "Neighbors considered must be above 1"
                    rlc::updateCharts(c("display_warning"))
                    return()
                  }

                  # Update how many neighbors to consider
                  n_NN <<- as.numeric(n_neighbors)

                  # Update k-scores and classification
                  avg_dist <<- Matrix::rowMeans(nn_subset$dist[,2:(n_NN+1)])

                  # Update display to be sure
                  display <<- "..."

                  # Update the corresponding charts
                  rlc::updateCharts()

                },
                place = "nn_considered",
                chartId = "input_nn_cons",
                width = 120,
                title = "Considered Neighbors",
                labels = c("")
  )

  # 5. Which gamma correction factor should be used (every chart will have to be updated)

  rlc::lc_input(type = "text",
                on_click = function(g) {

                  if (!(is.numeric(g))) {
                    display <<- "Enter a number"
                    rlc::updateCharts(c("display_warning"))
                    return()
                  }

                  # Update gamma
                  gamma <<- g

                  f_M_c <<- f_M^gamma

                  # Update the classification
                  classification <<- (f_M_c > threshold)

                  # Update charts
                  rlc::updateCharts()

                },
                rlc::dat(value = gamma),
                chartId = "input_gamma",
                place = "input_gamma",
                width = 120,
                title = "Gamma",
                labels = c("")
  )

  # 6. Which threshold to use to consider cells as marker expressing ones

  rlc::lc_input(type = "text",
                data = rlc::dat(value = round(threshold, 5)),
                on_click = function(thr) {

                  if (!(is.numeric(thr))) {
                    display <<- "Enter a number"
                    rlc::updateCharts(c("display_warning"))
                    return()
                  }

                  # Update the threshold
                  threshold <<- thr
                  display <<- "..."

                  # Update the classification
                  ifelse(above, classification <<- f_M_c > threshold, classification <<- f_M_c < threshold)

                  # Update charts
                  rlc::updateCharts()

                },
                chartId = "input_threshold",
                place = "input_threshold",
                width = 120,
                title = "Threshold",
                labels = c("")
  )

  # xxx. Change between above and below classification

  rlc::lc_input(type = "radio", labels = c("Above", "Below"),
                on_click = function(value) {

                  # Check which button was clicked and change to the corresponding embedding
                  if (value == 1) {

                    above <<- TRUE

                  } else if (value == 2) {

                    above <<- FALSE

                  }

                  # Classification
                  ifelse(above, classification <<- f_M_c > threshold, classification <<- f_M_c < threshold)

                  # Update the scatterplots
                  rlc::updateCharts()

                },
                value = 1,
                width = 120,
                place = "button_above_below",
                title = "Select"
  )

  # 11. Change between scatter plot and histogram to show fraction

  rlc::lc_input(type = "radio", labels = c("Scatter", "Histogram"),
                on_click = function(value) {

                  # Remove previous plot
                  rlc::removeChart("fraction_plot")

                  # Value == 1 if "Scatter" was clicked
                  if (value == 1) {

                    rlc::lc_scatter(
                      titleSize=16,
                      width=320,
                      height=320,
                      data = rlc::dat(
                        x = totals_ci^gamma,
                        y = f_M_c,
                        size = pt_size,
                        colourValue = classification,
                        colourDomain = c(TRUE, FALSE),
                        palette = c("green", "red")
                      ),
                      on_clickPosition = function(position) {

                        # Correct for the NA bug
                        if (!(is.na(position[1]))) {

                          #
                          x <- position[2]

                          # Update the threshold
                          threshold <<- x

                          # Update the classification
                          ifelse(above, classification <<- f_M_c > threshold, classification <<- f_M_c < threshold)

                          # Update charts
                          rlc::updateCharts()
                        }
                      },
                      transitionDuration = 0,
                      title = "Fraction vs. Library Size",
                      chartId = "fraction_plot",
                      place = "fraction_plot",
                      on_marked = function() {
                        rlc::mark(rlc::getMarked(chartId="fraction_plot", layerId = "Layer1"), "classification_plot", clear=TRUE)
                      }
                    )

                    rlc::lc_hLine(data = rlc::dat(
                      h = threshold
                    ),
                    addLayer = TRUE,
                    transitionDuration = 0,
                    chartId = "fraction_plot"
                    )

                  } else if (value == 2) {

                    # Histogram to show the fraction

                    rlc::lc_hist(
                      titleSize=16,
                      width=320,
                      height=320,
                      data = rlc::dat(
                        value = f_M_c,
                        nbins = 100
                      ),
                      on_clickPosition = function(position) {

                        # Get only the x position
                        x <- position[1]

                        # Correct for the NA bug
                        if (!(is.na(x))) {

                          # Update the threshold
                          threshold <<- x

                          # Update the classification
                          ifelse(above, classification <<- f_M_c > threshold, classification <<- f_M_c < threshold)

                          # Update charts
                          rlc::updateCharts()
                        }

                      },
                      title = "Fraction Histogram",
                      chartId = "fraction_plot",
                      transitionDuration = 0,
                      place = "fraction_plot"
                    )

                    rlc::lc_vLine(data = rlc::dat(
                      v = threshold
                    ),
                    addLayer = TRUE,
                    transitionDuration = 0,
                    chartId = "fraction_plot"
                    )

                  }

                },
                value = 1,
                width = 120,
                place = "button_hist_scatter_fraction",
                title = "Fraction Plot"
  )

  # 11. Change between different classification plots

  get_all_class_sample = function() {

    # Get cell types that were assigned in this sample
    types_sample <- names(g$rules[[sample_select]])

    # Create matrix to store the overall classification
    classification_matrix <- matrix(data=NA, nrow=sum(selection), ncol=1)
    rownames(classification_matrix) <- colnames(counts[,selection])
    colnames(classification_matrix) <- "classification"

    # Loop through all assigned cell types and transfer information to classification matrix
    for (type in types_sample) {
      # Get all the barcodes for the cell type when all == TRUE
      g$annotations[[sample_select]][[type]] %>% dplyr::filter(.data$all==TRUE) %>% dplyr::select(.data$barcode) %>% unlist -> barcodes_type
      # For now we don't know what the nature of the double assignment is
      classification_matrix[barcodes_type,1] <- ifelse(is.na(classification_matrix[barcodes_type,1]), type, "double")
    }
    return(classification_matrix %>% unlist)
  }

  # Some more helper functions, get the classification for the current cell type based on all applied rules
  helper_1 <- function() {
    return(g$annotations[[sample_select]][[current_class]]$all)
  }

  # Here is the actual input form

  rlc::lc_input(type = "radio", labels = c("Current Rule", "Current Cell Type", "All Cell Types"),
                on_click = function(value) {


                  rlc::removeChart("classification_plot")

                  if (value == 1) {

                    # 2D DimPlot to how show the classification based on rule

                    rlc::lc_scatter(
                      titleSize=16,
                      width=320,
                      height=320,
                      rlc::dat(
                        x = embedding[selection,1],
                        y = embedding[selection,2],
                        label = colnames(counts)[selection],
                        colourValue = classification,
                        colourDomain = c(TRUE, FALSE),
                        palette = c("green", "red"),
                        size = pt_size,
                        on_marked = function() {
                        }),
                      transitionDuration = 0,
                      title = "Classification Rule",
                      chartId = "classification_plot",
                      place = "classification_plot"
                    )

                  } else if (value == 2) {

                    # 2D DimPlot to how show the classification based on all rules for cell type

                    rlc::lc_scatter(
                      titleSize=16,
                      width=320,
                      height=320,
                      rlc::dat(
                        x = embedding[selection,1],
                        y = embedding[selection,2],
                        label = colnames(counts)[selection],
                        colourValue = helper_1(),
                        colourDomain = c(TRUE, FALSE),
                        palette = c("green", "red"),
                        size = pt_size,
                        on_marked = function() {
                        }),
                      transitionDuration = 0,
                      title = "Classification Cell Type",
                      chartId = "classification_plot",
                      place = "classification_plot"
                    )

                  } else if (value ==3) {

                    # 2D DimPlot to how show the classification based on all rules for cell type

                    rlc::lc_scatter(
                      titleSize=16,
                      width=320,
                      height=320,
                      rlc::dat(
                        x = embedding[selection,1],
                        y = embedding[selection,2],
                        label = colnames(counts)[selection],
                        colourValue = get_all_class_sample(),
                        #colourDomain = c(TRUE, FALSE),
                        #palette = c("green", "red"),
                        size = pt_size,
                        on_marked = function() {
                        }),
                      transitionDuration = 0,
                      title = "Classification Cell Type",
                      chartId = "classification_plot",
                      place = "classification_plot"
                    )

                  }

                },
                value = 1,
                width = 200,
                place = "button_classfification",
                title = "Classification Plot"
  )

  # 11. Change between different additional data plots

  rlc::lc_input(type = "radio", labels = c("Counts Marker", "Library Size"),
                on_click = function(value) {


                  rlc::removeChart("additional_data_plot")
                  # Check which button was clicked and change to the corresponding embedding
                  if (value == 1) {

                    # 8. 2D DimPlot to show the counts

                    rlc::lc_scatter(
                      titleSize=16,
                      width=320,
                      height=320,
                      rlc::dat(
                        x = embedding[selection,1],
                        y = embedding[selection,2],
                        label = colnames(counts)[selection],
                        colourValue = counts[marker,selection],
                        size = pt_size,
                        on_marked = function() {
                        }),
                      transitionDuration = 0,
                      title = "Counts",
                      chartId = "additional_data_plot",
                      place = "additional_data_plot"
                    )

                  } else if (value == 2) {

                    # 9. 2D DimPlot to show library size

                    rlc::lc_scatter(
                      titleSize=16,
                      width=320,
                      height=320,
                      rlc::dat(
                        x = embedding[selection,1],
                        y = embedding[selection,2],
                        label = colnames(counts)[selection],
                        colourValue = Matrix::colSums(counts[,selection]),
                        size = pt_size,
                        on_marked = function() {
                        }),
                      transitionDuration = 0,
                      title = "Library Size (total UMI count)",
                      chartId = "additional_data_plot",
                      place = "additional_data_plot"
                    )
                  }
                },
                value = 1,
                width = 200,
                place = "button_additional_data",
                title = "Additional Data Plot"
  )

  ###############################################################################

  # Define update functions

  # Update 1, which includes everything post sample selection
  # Will need to be called if the whole configuration is not already available
  update_1 = function() {

    # Get the counts for the marker genes for all cells within the selection.
    c_M <<- counts[marker, selection]
    counts_ci <<- sapply(1:nrow(nn_subset$idx), function(cell_i){
      neighbors <- nn_subset$idx[cell_i,1:k]
      sum(c_M[nn_subset$idx[cell_i,1:k]])

    })

    c_tot <<- Matrix::colSums(counts[, selection])
    totals_ci <<- sapply(1:nrow(nn_subset$idx), function(cell_i){
      sum(c_tot[nn_subset$idx[cell_i,1:k]])
    })

    # Calculate the fractions
    f_M <<- counts_ci/totals_ci
    f_M_c <<- f_M^gamma

    # Default positive classification if above threshold
    above <<- TRUE

    # Update the classification
    ifelse(above, classification <<- f_M_c > threshold, classification <<- f_M_c < threshold)

    # Update the display
    display <<- "..."

    # Update charts
    rlc::updateCharts()
  }


  # 7. Size of points in the scatter plots

  rlc::lc_input(type = "range",
                on_click = function(pt) {

                  # Update point size
                  pt_size <<- pt

                  # Update the pt size
                  rlc::updateCharts()

                },
                rlc::dat(value = pt_size),
                min = c(1),
                max = c(4),
                step=c(1),
                chartId = "input_pt",
                place = "slide_point_size",
                width = 240,
                title = "Point Size",
                labels = c("")
  )

  ###############################################################################

  # Add html output

  display <- "..."

  rlc::lc_html(rlc::dat(
    content = data.frame(Warning=display)),
    chartId = "display_warning",
    place = "warning"
  )

  # Include display to show classification
  current_class=NA   # initialize

  rlc::lc_html(rlc::dat(
    content = tibble::tibble(cell_type=current_class, gene=marker, above=above, gamma=gamma, threshold=threshold)
  ),
  place =  "display_rule"
  )

  # Include display to show all rules in sample

  get_all_rules = function(all_types=FALSE){
    cell_type <- names(g$rules[[sample_select]])

    display_tibble <- NULL

    for (name in cell_type) {
      display_tibble <- rbind(display_tibble,
                              g$rules[[sample_select]][[name]] %>%
                                dplyr::mutate(cell_type = name)) %>%
        dplyr::relocate(cell_type)
    }
    if (is.null(display_tibble)) {
      return("No Cell Type / Rule Created")
    } else if (all_types) {
      return(display_tibble)
    }
    else {
      return(display_tibble %>% dplyr::filter(cell_type==current_class) %>%
               dplyr::select(.data$gene, .data$gamma, .data$threshold, .data$above))
    }
  }

  rlc::lc_html(rlc::dat(
    content = get_all_rules(all_types = TRUE)
  ),
  place =  "rule_sample"
  )

  # Include display to show all rules for the current cell_type
  rlc::lc_html(rlc::dat(
    content = get_all_rules(all_types = FALSE)
  ),
  place =  "rule_celltype"
  )

  # Include display to show current cell

  rlc::lc_html(rlc::dat(
    content = tibble::tibble(Current_Cell = current_class)
  ),
  place =  "current_cell"
  )

  ###############################################################################

  # Create 2 list to store

  # a) annotations per cell type as tibble with barcodes, and rules as columns
  annotations <- vector(mode = "list", length=length(sample_info$sample))
  names(annotations) <- sample_info$sample

  # b) rules per cell type as names vector with gene, gamma, threshold, above
  rules <- vector(mode = "list", length=length(sample_info$sample))
  names(rules) <- sample_info$sample

  # Create empty list for each samples
  for (smp in sample_info$sample) {
    annotations[[smp]] <- vector(mode="list")
    rules[[smp]] <- vector(mode="list")
  }

  # Defining a new cell type
  rlc::lc_input(type = "text",
                width=240,
                data = rlc::dat(value = ifelse(is.na(current_class), "", current_class)),
                on_click = function(ct) {

                  if (!(is.character(ct))) {
                    display <<- "Enter a string"
                    rlc::updateCharts(c("display_warning"))
                    return()
                  } else if (ct == "") {
                    display <<- "Name cannot be empty"
                    rlc::updateCharts(c("display_warning"))
                    return()
                  }

                  current_class <<- ct

                },
                chartId = "input_class",
                place = "input_cell",
                title = "Define / Select Class",
                labels = c("")
  )

  # Adding a new cell

  rlc::lc_input(type = "button",
                width=100,
                labels = c("Add Cell"),
                on_click = function(value) {

                  rlc::updateCharts()

                },
                place = "button_cell")



  # Adding a new rule

  rlc::lc_input(type = "button",
           width=100,
           labels = c("Add Rule"),
           on_click = function(value) {

             # Only initiate a new tibble if the cell type has not been used
             if (!(current_class %in% names(g$annotations[[sample_select]]))) {

               # Initialize new tibble
               g$annotations[[sample_select]][[current_class]] <- tibble::tibble(
                 barcode = colnames(counts[,selection])
               )

               # Adding new column with the currently selected gene and its classification
               g$annotations[[sample_select]][[current_class]][marker] <- classification

               # Adding the meta data (what thresholds were used etc.)
               g$rules[[sample_select]][[current_class]] <- tibble::tibble(
                 gene = marker,
                 gamma = gamma,
                 threshold = threshold,
                 above = above
               )
             }

             # Check whether a rule has been created before for the given gene (in this case update and not append the tibble)
             if (marker %in% g$rules[[sample_select]][[current_class]]$gene) {
               g$rules[[sample_select]][[current_class]] <- g$rules[[sample_select]][[current_class]] %>%
                 dplyr::rows_update(tibble::tibble(
                   gene = marker,
                   gamma = gamma,
                   threshold = threshold,
                   above = above
                 ))

               # Append the new rule for the new gene
             } else
             {
               g$rules[[sample_select]][[current_class]] <- g$rules[[sample_select]][[current_class]] %>%
                 dplyr::add_row(
                   gene = marker,
                   gamma = gamma,
                   threshold = threshold,
                   above = above
                 )
             }

             g$annotations[[sample_select]][[current_class]]["all"] <- rep(NULL, length(classification))

             g$annotations[[sample_select]][[current_class]][marker] <- classification

             c_tmp <- colnames(g$annotations[[sample_select]][[current_class]])

             g$annotations[[sample_select]][[current_class]] <- g$annotations[[sample_select]][[current_class]] %>%
               dplyr::rowwise() %>%
               # For checking all, take all columns but barcodes and all column (found thus far no easier way!)
               dplyr::mutate(all = all(dplyr::c_across(cols = setdiff(c_tmp, c("barcode", "all")))))
             # There must be simpler ways though!

             rlc::updateCharts()

           },
           place = "button_rule")



  # Delete rule

  rlc::lc_input(type = "text",
                width=110,
                labels = c(""),
                title = "Delete Rule for Gene:",
                on_click = function(value) {

                  if(value %in% g$rules[[sample_select]][[current_class]]["gene"]) {

                    g$rules[[sample_select]][[current_class]] <- g$rules[[sample_select]][[current_class]] %>%
                      dplyr::filter(.data$gene != value)

                    g$annotations[[sample_select]][[current_class]] <- g$annotations[[sample_select]][[current_class]] %>%
                      dplyr::select(-dplyr::matches(value))

                    # Check whether rules are left
                    if (length(setdiff(colnames(g$annotations[[sample_select]][[current_class]]), c("barcode", "all"))) == 0) {
                      g$annotations[[sample_select]][[current_class]] <- g$annotations[[sample_select]][[current_class]] %>%
                        dplyr::mutate(all = FALSE)
                    } else {
                      # Reevaluate all columns
                      g$annotations[[sample_select]][[current_class]] <- g$annotations[[sample_select]][[current_class]] %>%
                        dplyr::rowwise() %>%
                        # Check all column but "barcode" and "all", wether the entries in the row are True
                        dplyr::mutate(all = all(dplyr::c_across(cols = setdiff(colnames(.data), c("barcode", "all")))))
                    }

                    # Update everything
                    update_1()
                  } else {
                    display <<- "No rule for this gene"
                  }
                  rlc::updateCharts()
                },
                place = "delete_rule")
}
