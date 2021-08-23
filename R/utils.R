#' @importFrom dplyr %>%
#' @export
dplyr::`%>%`

#' Generate sample info data frame
#' @name generate_sample_info
#' @importFrom dplyr %>%
#' @export
generate_sample_info <- function(cells) {
  tibble::tibble(sample = unique(cells$sample)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(total = sum(cells$sample == sample)) %>%
    dplyr::ungroup() -> sample_info
  
  # Correction factor adjust indices for nn below
  sample_info$correction = c(0, cumsum(sample_info$total)[1:nrow(sample_info)-1])
  return(sample_info)
}

#'
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
  return(list("marker_counts" = counts_ci,
           "total_counts" = totals_ci,
           "smoothed_fraction" = f_M,
           "corrected_fraction" = f_M_c))
}

#' 
#' @name smoothed_fraction_wrapper
smoothed_fraction_wrapper = function(app_env) {
  tmp_calc <- smoothed_fraction(app_env$counts, app_env$marker, app_env$selection, app_env$nn_subset, app_env$gamma)
  app_env$counts_ci <- tmp_calc$marker_counts
  app_env$totals_ci <- tmp_calc$total_counts
  app_env$f_M <- tmp_calc$smoothed_fraction
  app_env$f_M_c <- tmp_calc$corrected_fraction
}

#' Generate nearest neighbor per sample
#' @name run_nn
#' @importFrom dplyr %>%
#' @param cells data.frame with cell metadata
#' @param pc_space PCA coordinates
#' @param k number of neighbors considered for smoothing
#' @param dim number of PC dimensions used for nearest neighbor graph
#' @export
run_nn <- function(cells, pc_space, k=50, dim=30) {

  tibble::tibble(sample = unique(cells$sample)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(total = sum(cells$sample == sample)) %>%
    dplyr::ungroup() -> sample_info

  # Correction factor adjust indices for nn below
  sample_info$correction = c(0, cumsum(sample_info$total)[1:nrow(sample_info)-1])

  # tmp variables
  idx_tmp <- NULL
  dist_tmp <- NULL

  for (smp in sample_info$sample) {
    ifelse(ncol(pc_space) < dim,
           nn_tmp <- RANN::nn2(pc_space[cells$sample==smp,], k=k),
           nn_tmp <- RANN::nn2(pc_space[cells$sample==smp,1:dim], k=k)
           )
    correction <- sample_info %>% dplyr::filter(sample==smp) %>% dplyr::pull(correction)
    idx_tmp <- rbind(idx_tmp, nn_tmp$nn.idx + correction)
    dist_tmp <- rbind(dist_tmp, nn_tmp$nn.dists)
  }
  nn <- list(
    idx = idx_tmp,
    dists = dist_tmp
  )
  stopifnot(nn$idx[,1] == 1:nrow(cells))
  return(nn)
}


#' Generate nearest neighbor per sample
#' @name run_classification
#' @importFrom dplyr %>%
#' @param annotations list created by lc_vis
run_classification = function(annotations, sample_select) {
  # Get cell types that were assigned in this sample
  types_sample <- names(annotations[[sample_select]])
  
  # Initialize tibble to store classifications
  tmp_tibble <- tibble::tibble(
    barcode = annotations[[sample_select]][[1]]$barcode,
    classification = "none")
  
  for (type in types_sample) {
    
    tmp_tibble <- tmp_tibble %>%
      dplyr::mutate(tmp = annotations[[sample_select]][[type]]$all) %>%
      dplyr::mutate(classification = dplyr::case_when(
        # ambigious if classification is not nonce and another class is positive for the cell
        !(classification=="none") & .data$tmp ~ "ambigious",
        # if classification is nonce and the cell is positive for the class assign the type
        classification=="none" & .data$tmp ~ type,
        TRUE ~ classification
      )) }
  
  return(tmp_tibble %>% dplyr::pull(classification))
}

#' Generate nearest neighbor per sample
#' @name run_classification_wrapper
#' @importFrom dplyr %>%
#' @param annotations list created by lc_vis
#' @export
run_classification_wrapper = function(annotations) {

  annotated_samples <- names(annotations)
  ann_list <- vector(mode = "list", length=length(annotated_samples))

  for (smp in annotated_samples) {

    # Check whether there are actually annotations
    if (length(annotations[[smp]])>0) {
      # Get cell types that were assigned in this sample
      types_sample <- names(annotations[[smp]])

      tmp_tibble <- tibble::tibble(
        barcode = annotations[[smp]][[1]]$barcode,
        classification = "none")

      for (type in types_sample) {

        # type_annotation should be a logical vector
        type_annotation <- annotations[[smp]][[type]]$all

        tmp_tibble <- tmp_tibble %>%
          dplyr::mutate(tmp = type_annotation) %>%
          dplyr::mutate(classification = dplyr::case_when(
            # ambigious if classification is not NA and another class is positive for the cell
            !(classification=="none") & .data$tmp ~ "ambigious",
            # if classification is NA and the cell is positive for the class assign the type
            classification=="none" & .data$tmp ~ type,
            TRUE ~ classification
          ))

        ann_list[[smp]] <- tmp_tibble %>% dplyr::select(-.data$tmp)
      }
    }
  }
  return(ann_list)
}

#' 
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @name get_all_rules
#' @export
get_all_rules = function(rules, sample_select, all_types, current_class){
  cell_type <- names(rules[[sample_select]])
  
  display_tibble <- NULL
  
  for (name in cell_type) {
    display_tibble <- rbind(display_tibble,
                            rules[[sample_select]][[name]] %>%
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

#' 
#' 
initialize_app_env = function(app_env, cells, counts, pc_space, embedding, nn, k) {
  app_env$cells <- cells
  app_env$counts <- counts
  app_env$pc_space <- pc_space
  app_env$embedding <- embedding
  app_env$nn <- nn
}