#' Update classification
#' @name update_classification
update_classification = function(app_env) {
  ifelse(app_env$above, 
         app_env$classification <- app_env$f_M_c > app_env$threshold, 
         app_env$classification <- app_env$f_M_c < app_env$threshold)

  # Update warning display and all charts
  app_env$display <- "..."
}

#' function to calculate the fraction smoothed over nearest neighbors
#' @name smoothed_fraction
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
smoothed_fraction_wrapper = function(app_env) {
  tmp_calc <- smoothed_fraction(app_env$counts, app_env$marker, app_env$selection, app_env$nn_subset, app_env$gamma)
  app_env$counts_ci <- tmp_calc$marker_counts
  app_env$totals_ci <- tmp_calc$total_counts
  app_env$f_M <- tmp_calc$smoothed_fraction
  app_env$f_M_c <- tmp_calc$corrected_fraction
}

#' Update upon selection of new sample
#' @name update_new_sample
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