#' Update classification
#' @name update_classification
#' 
#' @param app_env R environment, app environment
#' 
update_classification = function(app_env) {
  ifelse(app_env$above, 
         app_env$classification <- app_env$f_M_c > app_env$threshold, 
         app_env$classification <- app_env$f_M_c < app_env$threshold)

  # Update warning display and all charts
  app_env$display <- "..."
}

#' function to calculate the fraction smoothed over nearest neighbors
#' @name smoothed_fraction
#' 
#' @param counts sparse matrix, containing counts (cols = cells, genes = rows)
#' @param marker string, selected marker gene
#' @param selection logical array, used to subset cells belonging to the selected sample
#' @param nn_subset list, nearest neighbor indices and distances for cells in the sample
#' @param gamma numeric, gamma-correction factor
#' 
smoothed_fraction <- function(counts, marker, selection, nn_subset, gamma) {
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
  f_M_c <- f_M^gamma
  
  # Return everything in one list
  return(list("marker_counts" = counts_ci,
              "total_counts" = totals_ci,
              "smoothed_fraction" = f_M,
              "corrected_fraction" = f_M_c))
}

#' Wrapper calling smoothed_fraction function and adding the result to the app environment
#' @name smoothed_fraction_wrapper
#' 
#' @param app_env R environment, app enviroment
#' 
smoothed_fraction_wrapper = function(app_env) {
  tmp_calc <- smoothed_fraction(app_env$counts, app_env$marker, app_env$selection, app_env$nn_subset, app_env$gamma)
  app_env$counts_ci <- tmp_calc$marker_counts
  app_env$totals_ci <- tmp_calc$total_counts
  app_env$f_M <- tmp_calc$smoothed_fraction
  app_env$f_M_c <- tmp_calc$corrected_fraction
}

#' Update upon selection of new sample
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @name update_new_sample
#' 
#' @param app_env R environment, app enviroment
#' @param smp string, selected sample
#' 
update_new_sample = function(app_env, smp) {
  
  # Update selected sample and subset
  app_env$sample_select <- smp
  app_env$selection <- app_env$cells$sample==app_env$sample_select
  
  # Getting the nearest neighbors of all cells in the subset.
  app_env$correction <- app_env$sample_info %>% dplyr::filter(sample==app_env$sample_select) %>% dplyr::pull(.data$correction)
  app_env$nn_subset <- list(idx = (app_env$nn$idx[app_env$selection,] - app_env$correction),
                     dists = app_env$nn$dist[app_env$selection,])
  
  # Get the average distance for cell_i to x nearest neighbors
  app_env$avg_dist <- Matrix::rowMeans(app_env$nn_subset$dist[,2:(app_env$n_NN+1)])
}

#' Function called if new class is created in a sample
#' @name udpate_new_class
#' 
#' @param app_env R environment, app enviroment
#' @param g_env R environment, global environment
#' 
update_new_class <- function(g_env, app_env) {
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

#' Function called if rule for gene exists and thus needs to be updated
#' @importFrom dplyr %>%
#'
#' @name update_rule
#'
#' @param app_env R environment, app enviroment
#' @param g_env R environment, global environment
#'
update_rule <- function(g_env, app_env) {
  g_env$rules[[app_env$sample_select]][[app_env$current_class]] <- g_env$rules[[app_env$sample_select]][[app_env$current_class]] %>%
    dplyr::rows_update(tibble::tibble(
      gene = app_env$marker,
      gamma = app_env$gamma,
      threshold = app_env$threshold,
      above = app_env$above
    ))
}

#' Function called if rule for gene did not exist and thus needs to be added
#' @importFrom dplyr %>%
#'
#' @name add_rules
#'
#' @param app_env R environment, app enviroment
#' @param g_env R environment, global environment
#'
add_rules <- function(g_env, app_env) {
  g_env$rules[[app_env$sample_select]][[app_env$current_class]] <- g_env$rules[[app_env$sample_select]][[app_env$current_class]] %>%
    dplyr::add_row(
      gene = app_env$marker,
      gamma = app_env$gamma,
      threshold = app_env$threshold,
      above = app_env$above
    )
}

#' Function called after rules was added/deleted to update the cell type annotation
#' @importFrom dplyr %>%
#'
#' @name update_annotation
#'
#' @param app_env R environment, app enviroment
#' @param g_env R environment, global environment
#'
update_annotation <- function(g_env, app_env) {
  
  c_names <- colnames(g_env$annotations[[app_env$sample_select]][[app_env$current_class]])
  
  genes <- setdiff(c_names, c("barcode", "all"))
  
  if (length(genes) == 0) {
    g_env$annotations[[app_env$sample_select]][[app_env$current_class]]["all"] <- FALSE
  } else {
    # For checking whether a cell adheres to all rules we have to check whether all column (=genes) are TRUE
    tmp <- g_env$annotations[[app_env$sample_select]][[app_env$current_class]] %>%
      dplyr::select(genes)
    
    all_classification <- apply(tmp, MARGIN=1, FUN=all)
    
    g_env$annotations[[app_env$sample_select]][[app_env$current_class]]["all"] <- all_classification
  }
}





