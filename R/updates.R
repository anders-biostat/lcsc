#' Update classification
#' @name update_classification
update_classification = function(app_env) {
  ifelse(app_env$above, 
         app_env$classification <- app_env$f_M_c > app_env$threshold, 
         app_env$classification <- app_env$f_M_c < app_env$threshold)

  # Update warning display and all charts
  app_env$display <- "..."
}

#' Update upon selection of new sample
#' @name update_new_sample
update_new_sample = function(app_env, smp) {
  
  # Update selected sample and subset
  app_env$sample_select <- smp
  app_env$selection <- app_env$cells$sample==app_env$sample_select
  
  # Getting the nearest neighbors of all cells in the subset.
  app_env$correction <- app_env$sample_info %>% dplyr::filter(sample==app_env$sample_select) %>% dplyr::pull(correction)
  app_env$nn_subset <- list(idx = (app_env$nn$idx[app_env$selection,] - app_env$correction),
                     dists = app_env$nn$dist[app_env$selection,])
  
  # Get the average distance for cell_i to x nearest neighbors
  app_env$avg_dist <- Matrix::rowMeans(app_env$nn_subset$dist[,2:(app_env$n_NN+1)])
}