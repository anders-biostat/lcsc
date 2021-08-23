
#' Visualize and annotate single cell data
#' @name lc_vis
#' 
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' 
#' @param cells data.frame with cell metadata
#' @param counts sparse count matrix (cols = cells, genes = rows)
#' @param pc_space PCA coordinates (cols = coordiantes, rows = cells)
#' @param embedding coordinates for two-dimensional embedding (cols = coordiantes, rows = cells)
#' @param nn list containing nearest neighbor information (created by run_nn function)
#' @param k number of neighbors considered for smoothing (default k = 50)
#' 
#' @export
lc_vis <- function(cells, counts, pc_space, embedding, nn, k=50){

  ###############################################################################
  
  g_env <- globalenv()
  
  # Generate and populate app environment
  app_env <- new.env()
  initialize_app_env(app_env, cells, counts, pc_space, embedding, nn, k)
  rm(cells, counts, pc_space, embedding, nn, k)
  
  app_env$sample_info <- generate_sample_info(app_env$cells)

  ###############################################################################

  # Start rlc page
  startPage = system.file("http_root/rlc_layout.html", package="lcsc")
  rlc::openPage(startPage=startPage, useViewer=FALSE)

  ###############################################################################

  # Initialize the application choosing a random gene and the first sample
  app_env$marker <- sample(rownames(app_env$counts),1)
  app_env$sample_select <- app_env$sample_info$sample[1]
  app_env$selection <- cells$sample==app_env$sample_select

  # Getting the nearest neighbors of all cells in the subset.
  app_env$correction <- app_env$sample_info %>% dplyr::filter(sample==app_env$sample_select) %>% dplyr::pull(correction)
  app_env$nn_subset <- list(idx = (app_env$nn$idx[app_env$selection,] - app_env$correction),
                    dists = app_env$nn$dist[app_env$selection,])

  # Default parameters
  app_env$n_NN <- 25
  app_env$gamma = 0.5
  app_env$threshold <- 0.02
  app_env$above <- TRUE
  
  # Get the average distance for cell_i to n_NN nearest neighbors.
  app_env$avg_dist <- Matrix::rowMeans(app_env$nn_subset$dists[,2:(app_env$n_NN+1)])

  # Calculate the marker gene fraction smoothed over the neighbors
  smoothed_fraction_wrapper(app_env)
  
  # Classify cells according to whether the smoothed fraction is above or below the threshold
  update_classification(app_env)

  ###############################################################################

  # Initial constant for figures
  app_env$pt_size <- 1 # point size

  ###############################################################################

  # Initialize the plots

  # 1. Histogram to get feeling for the distances between the cells

  rlc::lc_hist(
    titleSize=16,
    width=320,
    height=320,
    data = rlc::dat(
      value = app_env$avg_dist,
      nbins = 100,
      title = paste0("Mean dist to ", app_env$n_NN, "th nearest neigbor")
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
      x = app_env$embedding[app_env$selection,1],
      y = app_env$embedding[app_env$selection,2],
      label = colnames(app_env$counts)[app_env$selection],
      colourValue = app_env$avg_dist,
      size = app_env$pt_size,
      title = paste0("Mean dist to ", app_env$n_NN, "th nearest neigbor"),
      on_marked = function() {
      }),
    transitionDuration = 0,
    chartId = "nn_below_plot",
    place = "nn_below_plot"
  )

  # 3. Plots to show the distribution of the marker gene fraction

  rlc::lc_scatter(
    titleSize=16,
    width=320,
    height=320,
    data = rlc::dat(
      x = app_env$totals_ci^app_env$gamma,
      y = app_env$f_M_c,
      size = app_env$pt_size,
      colourValue = app_env$classification,
      colourDomain = c(TRUE, FALSE),
      palette = c("green", "red")
    ),
    on_clickPosition = function(position) {

      # Correct for the NA bug (on_ClickPosition return first NA and then the value)
      if (!(is.na(position[1]))) {

        # Update the threshold, the classification and the charts
        app_env$threshold <- position[2]
        update_classification(app_env)
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
    h = app_env$threshold
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
      x = app_env$embedding[app_env$selection,1],
      y = app_env$embedding[app_env$selection,2],
      label = colnames(app_env$counts)[app_env$selection],
      colourValue = app_env$f_M_c,
      size = app_env$pt_size
      ),
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
      x = app_env$embedding[app_env$selection,1],
      y = app_env$embedding[app_env$selection,2],
      label = colnames(app_env$counts)[app_env$selection],
      colourValue = app_env$classification,
      colourDomain = c(TRUE, FALSE),
      palette = c("green", "red"),
      size = app_env$pt_size,
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
      x = app_env$embedding[app_env$selection,1],
      y = app_env$embedding[app_env$selection,2],
      label = colnames(app_env$counts)[app_env$selection],
      colourValue = app_env$counts[app_env$marker,app_env$selection],
      size = app_env$pt_size),
    transitionDuration = 0,
    title = "Counts",
    chartId = "additional_data_plot",
    place = "additional_data_plot"
  )

  ###############################################################################

  # Define the input

  # 1. Which sample is used.

  rlc::lc_input(type = "text", label=c(""),
                rlc::dat(value = app_env$sample_select),
                on_click = function(smp) {

                  # Add security check, whether the gene really exists
                  if (!(smp %in% app_env$sample_info$sample)) {
                    
                    app_env$display <- "Sample does not exist"
                    rlc::updateChartsc(c("display_warning"))
                    
                  } else {
                  
                  update_new_sample(app_env, smp)

                  smoothed_fraction_wrapper(app_env)
                  
                  update_classification(app_env)
                  
                  rlc::updateCharts() }

                },
                place = "sample",
                chartId = "input_sample",
                width = 240,
                title = "Sample"
  )

  # 2. Which marker gene is used.

  rlc::lc_input(type = "text", label=c(""),
                rlc::dat(value = app_env$marker),
                on_click = function(gene) {

                  # Add security check, whether the gene really exists
                  if (!(gene %in% rownames(app_env$counts))) {
                    
                    app_env$display <- "Gene does not exist"
                    rlc::updateCharts(c("display_warning"))
                    
                  } else {

                  # Update marker
                  app_env$marker <- gene

                  smoothed_fraction_wrapper(app_env)
                  
                  update_classification(app_env) }
                  
                  rlc::updateCharts()

                },
                place = "gene",
                chartId = "input_marker",
                width = 240,
                title = "Marker Gene"
  )

  # 3. How many neighbors to consider to show average distance

  rlc::lc_input(type = "text",
                rlc::dat(value = app_env$n_NN),
                on_click = function(n_neighbors) {

                  if (!(is.numeric(n_neighbors))) {
                    app_env$display <- "Enter a number"
                    rlc::updateCharts(c("display_warning"))
                    return()
                  }

                  # Checking
                  if (n_neighbors >= k) {
                    app_env$display <- paste("Neighbors considered must be below", k)
                    rlc::updateCharts(c("display_warning"))
                    return()
                  } else if (n_neighbors < 2) {
                    app_env$display <- "Neighbors considered must be above 1"
                    rlc::updateCharts(c("display_warning"))
                    return()
                  }

                  # Update how many neighbors to consider
                  app_env$n_NN <- as.numeric(n_neighbors)
                  app_env$avg_dist <- Matrix::rowMeans(app_env$nn_subset$dist[,2:(app_env$n_NN+1)])
                  app_env$display <- "..."

                  # Update the corresponding charts
                  rlc::updateCharts()

                },
                place = "nn_considered",
                chartId = "input_nn_cons",
                width = 120,
                title = "Considered Neighbors",
                labels = c("")
  )

  # 4. Which gamma correction factor should be used (every chart will have to be updated)

  rlc::lc_input(type = "text",
                on_click = function(g) {

                  if (!(is.numeric(g))) {
                    app_env$display <- "Enter a number"
                    rlc::updateCharts(c("display_warning"))
                    return()
                  }

                  # Update gamma
                  app_env$gamma <- g
                  app_env$f_M_c <- app_env$f_M^app_env$gamma

                  # Update the classification
                  update_classification(app_env)

                  # Update charts
                  rlc::updateCharts()

                },
                rlc::dat(value = app_env$gamma),
                chartId = "input_gamma",
                place = "input_gamma",
                width = 120,
                title = "Gamma",
                labels = c("")
  )

  # 5. Which threshold to use to consider cells as marker expressing ones

  rlc::lc_input(type = "text",
                data = rlc::dat(value = round(app_env$threshold, 5)),
                on_click = function(thr) {

                  if (!(is.numeric(thr))) {
                    app_env$display <- "Enter a number"
                    rlc::updateCharts(c("display_warning"))
                    return()
                  }

                  # Update the threshold and the 
                  app_env$threshold <- thr
                  update_classification(app_env)
                  
                  # Update charts
                  rlc::updateCharts()

                },
                chartId = "input_threshold",
                place = "input_threshold",
                width = 120,
                title = "Threshold",
                labels = c("")
  )

  # 6. Change between above and below classification

  rlc::lc_input(type = "radio", labels = c("Above", "Below"),
                on_click = function(value) {

                  # Check which button was clicked and change to the corresponding embedding
                  if (value == 1) { app_env$above <- TRUE
                  } else if (value == 2) { app_env$above <- FALSE
                  }

                  update_classification(app_env)
                  
                  # Update the scatterplots
                  rlc::updateCharts()

                },
                value = 1,
                width = 120,
                place = "button_above_below",
                title = "Select"
  )

  # 7. Change between scatter plot and histogram to show fraction

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
                        x = app_env$totals_ci^app_env$gamma,
                        y = app_env$f_M_c,
                        size = app_env$pt_size,
                        colourValue = app_env$classification,
                        colourDomain = c(TRUE, FALSE),
                        palette = c("green", "red")
                      ),
                      on_clickPosition = function(position) {

                        # Correct for the rlc NA bug
                        if (!(is.na(position[1]))) {

                          # Update the threshold
                          app_env$threshold <- position[2]

                          update_classification(app_env)
                          
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
                      h = app_env$threshold
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
                        value = app_env$f_M_c,
                        nbins = 100
                      ),
                      on_clickPosition = function(position) {

                        # Correct for the rlc NA bug
                        if (!(is.na(position[1]))) {

                          # Update the threshold
                          app_env$threshold <- position[1]

                          update_classification(app_env)
                          
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
                      v = app_env$threshold
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

  # 8. Change between different classification plots

  rlc::lc_input(type = "radio", 
                labels = c("Current Rule", "Current Cell Type", "All Cell Types"),
                on_click = function(value) {

                  rlc::removeChart("classification_plot")

                  if (value == 1) {

                    # 2D DimPlot to how show the classification based on rule

                    rlc::lc_scatter(
                      titleSize=16,
                      width=320,
                      height=320,
                      rlc::dat(
                        x = app_env$embedding[app_env$selection,1],
                        y = app_env$embedding[app_env$selection,2],
                        label = colnames(app_env$counts)[app_env$selection],
                        colourValue = app_env$classification,
                        colourDomain = c(TRUE, FALSE),
                        palette = c("green", "red"),
                        size = app_env$pt_size,
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
                        x = app_env$embedding[app_env$selection,1],
                        y = app_env$embedding[app_env$selection,2],
                        label = colnames(app_env$counts)[app_env$selection],
                        colourValue = g_env$annotations[[app_env$sample_select]][[app_env$current_class]]$all,
                        colourDomain = c(TRUE, FALSE),
                        palette = c("green", "red"),
                        size = app_env$pt_size,
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
                        x = app_env$embedding[app_env$selection,1],
                        y = app_env$embedding[app_env$selection,2],
                        label = colnames(app_env$counts)[app_env$selection],
                        colourValue = run_classification(g_env$annotations, app_env$sample_select),
                        size = app_env$pt_size,
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

  # 9. Change between different additional data plots

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
                        x = app_env$embedding[app_env$selection,1],
                        y = app_env$embedding[app_env$selection,2],
                        label = colnames(app_env$counts)[app_env$selection],
                        colourValue = app_env$counts[app_env$marker,app_env$selection],
                        size = app_env$pt_size),
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
                        x = app_env$embedding[app_env$selection,1],
                        y = app_env$embedding[app_env$selection,2],
                        label = colnames(app_env$counts)[app_env$selection],
                        colourValue = Matrix::colSums(app_env$counts[,app_env$selection]),
                        size = app_env$pt_size),
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

  # 10. Size of points in the scatter plots

  rlc::lc_input(type = "range",
                on_click = function(pt) {

                  # Update point size
                  app_env$pt_size <- pt
                  rlc::updateCharts()

                },
                rlc::dat(value = app_env$pt_size),
                min = c(1),
                max = c(4),
                step=c(0.5),
                chartId = "input_pt",
                place = "slide_point_size",
                width = 240,
                title = "Point Size",
                labels = c("")
  )

  ###############################################################################

  # Add html output

  app_env$display <- "..."

  # 1. Display to show warnings
  
  rlc::lc_html(rlc::dat(
    content = data.frame(Warning=app_env$display)),
    chartId = "display_warning",
    place = "warning"
  )

  # 2. Display to show classification
  
  app_env$current_class=NA   # initialize

  rlc::lc_html(rlc::dat(
    content = tibble::tibble(cell_type=app_env$current_class, 
                             gene=app_env$marker, above=app_env$above, 
                             gamma=app_env$gamma, threshold=app_env$threshold)
  ),
  place =  "display_rule"
  )

  # 3. Display to show all rules in sample
  
  rlc::lc_html(rlc::dat(
    content = get_all_rules(g_env$rules, app_env$sample_select, all_types=TRUE, app_env$current_class)
  ),
  place =  "rule_sample"
  )

  # 4. Display to show all rules for the current cell_type
  
  rlc::lc_html(rlc::dat(
    content = get_all_rules(g_env$rules, app_env$sample_select, all_types=FALSE, app_env$current_class)
  ),
  place =  "rule_celltype"
  )

  # 5. Display to show current cell

  rlc::lc_html(rlc::dat(
    content = tibble::tibble(Current_Cell = app_env$current_class)
  ),
  place =  "current_cell"
  )

  ###############################################################################

  # Create 2 list to store

  # a) annotations per cell type as tibble with barcodes, and rules as columns
  g_env$annotations <- vector(mode = "list", length=length(app_env$sample_info$sample))
  names(g_env$annotations) <- app_env$sample_info$sample

  # b) rules per cell type as names vector with gene, gamma, threshold, above
  g_env$rules <- vector(mode = "list", length=length(app_env$sample_info$sample))
  names(g_env$rules) <- app_env$sample_info$sample

  # Create empty list for each samples
  for (smp in app_env$sample_info$sample) {
    g_env$annotations[[smp]] <- vector(mode="list")
    g_env$rules[[smp]] <- vector(mode="list")
  }

  # Defining a new cell type
  rlc::lc_input(type = "text",
                width=240,
                data = rlc::dat(
                  value = ifelse(is.na(app_env$current_class), "", app_env$current_class)),
                on_click = function(ct) {

                  if (!(is.character(ct))) {
                    app_env$display <- "Enter a string"
                    rlc::updateCharts(c("display_warning"))
                    return()
                  } else if (ct == "") {
                    app_env$display <- "Name cannot be empty"
                    rlc::updateCharts(c("display_warning"))
                    return()
                  }

                  app_env$current_class <- ct

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

             if (is.na(app_env$current_class)) {
               app_env$display <- "Cannot add rule without class"
               rlc::updateCharts(c("display_warning"))
               return()
             } else {
               app_env$display <- "..."
               rlc::updateCharts(c("display_warning"))
             }

             # Only initiate a new tibble if the cell type has not been used
             if (!(app_env$current_class %in% names(g_env$annotations[[app_env$sample_select]]))) {

               # Initialize new tibble
               g_env$annotations[[app_env$sample_select]][[app_env$current_class]] <- tibble::tibble(
                 barcode = colnames(app_env$counts[,app_env$selection])
               )

               # Adding new column with the currently selected gene and its classification
               g_env$annotations[[app_env$sample_select]][[app_env$current_class]][app_env$marker] <- app_env$classification

               # Adding the rule data
               g_env$rules[[app_env$sample_select]][[app_env$current_class]] <- tibble::tibble(
                 gene = app_env$marker,
                 gamma = app_env$gamma,
                 threshold = app_env$threshold,
                 above = app_env$above
               )
             }

             # Check whether a rule has been created before for the given gene
             if (app_env$marker %in% g_env$rules[[app_env$sample_select]][[app_env$current_class]]$gene) # update the tibble
               {
               g_env$rules[[app_env$sample_select]][[app_env$current_class]] <- g_env$rules[[app_env$sample_select]][[app_env$current_class]] %>%
                 dplyr::rows_update(tibble::tibble(
                   gene = app_env$marker,
                   gamma = app_env$gamma,
                   threshold = app_env$threshold,
                   above = app_env$above
                 ))


             } else # Append the new rule for the new gene
             {
               g_env$rules[[app_env$sample_select]][[app_env$current_class]] <- g_env$rules[[app_env$sample_select]][[app_env$current_class]] %>%
                 dplyr::add_row(
                   gene = app_env$marker,
                   gamma = app_env$gamma,
                   threshold = app_env$threshold,
                   above = app_env$above
                 )
             }

             # Initializing or resetting the all column
             g_env$annotations[[app_env$sample_select]][[app_env$current_class]]["all"] <- rep(NULL, length(app_env$classification))

             g_env$annotations[[app_env$sample_select]][[app_env$current_class]][app_env$marker] <- app_env$classification

             c_names <- colnames(g_env$annotations[[app_env$sample_select]][[app_env$current_class]])

             # For checking whether a cell adheres to all rules we have to check all column (=genes) are TRUE
             g_env$annotations[[app_env$sample_select]][[app_env$current_class]] <- g_env$annotations[[app_env$sample_select]][[app_env$current_class]] %>%
               dplyr::rowwise() %>%
               dplyr::mutate(all = all(dplyr::c_across(cols = setdiff(c_names, c("barcode", "all"))))) # select all columns but barcode & all

             rlc::updateCharts()

           },
           place = "button_rule")



  # Delete rule

  rlc::lc_input(type = "text",
                width=110,
                labels = c(""),
                title = "Delete Rule for Gene:",
                on_click = function(value) {

                  if(value %in% dplyr::pull(g_env$rules[[app_env$sample_select]][[app_env$current_class]], .data$gene)) {

                    g_env$rules[[app_env$sample_select]][[app_env$current_class]] <- g_env$rules[[app_env$sample_select]][[app_env$current_class]] %>%
                      dplyr::filter(.data$gene != value)

                    g_env$annotations[[app_env$sample_select]][[app_env$current_class]] <- g_env$annotations[[app_env$sample_select]][[app_env$current_class]] %>%
                      dplyr::select(-dplyr::matches(value))

                    # Check whether rules are left
                    if (length(setdiff(colnames(g_env$annotations[[app_env$sample_select]][[app_env$current_class]]), c("barcode", "all"))) == 0) {
                      g_env$annotations[[app_env$sample_select]][[app_env$current_class]] <- g_env$annotations[[app_env$sample_select]][[app_env$current_class]] %>%
                        dplyr::mutate(all = FALSE)
                    } else # Reevaluate all columns as done above
                      {
                      g_env$annotations[[app_env$sample_select]][[app_env$current_class]] <- g_env$annotations[[app_env$sample_select]][[app_env$current_class]] %>%
                        dplyr::rowwise() %>%
                        dplyr::mutate(all = all(dplyr::c_across(cols = setdiff(colnames(.data), c("barcode", "all")))))
                    }

                    # Update
                    smoothed_fraction_wrapper(app_env)
                    update_classification(app_env)
                    rlc::updateCharts()
                    
                    
                  } else {
                    print("evaluates to FALSE")
                    app_env$display <- "No rule for this gene"
                  }
                  rlc::updateCharts()
                },
                place = "delete_rule")
  
  # tmp
  g_env$tmp <- app_env
}


