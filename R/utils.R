#' Generate sample info data frame
#' @name generate_sample_info
#' @importFrom dplyr %>%
#' @param cells tibble, containing barcodes and corresponding sample
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

#' Generate nearest neighbor per sample
#' @name run_nn
#' 
#' @importFrom dplyr %>%
#' 
#' @param cells data.frame, containing cell metadata
#' @param pc_space matrix, containing PCA coordinates (cols = coordiantes, rows = cells)
#' @param k numeric, number of neighbors considered for smoothing (default k = 50)
#' @param dim, numeric, number of PC dimensions used to compute nearest neighbor graph
#' 
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
#' 
#' @importFrom dplyr %>%
#' 
#' @param annotations list, created by lc_vis
#' @param sample_select string, selected sample which must be in sample_info$sample
#' @param app_env R environment, app environment
#' 
run_classification = function(annotations, sample_select, app_env) {
  
  # Initialize tibble to store classifications
  tmp_tibble <- tibble::tibble(
    barcode = colnames(app_env$counts[,app_env$selection]),
    classification = "none")
    
  # Get cell types that were assigned in this sample
  types_sample <- names(annotations[[sample_select]])
  
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
  
  return(tmp_tibble %>% dplyr::pull(.data$classification))
}

#' Generate nearest neighbor per sample
#' @name classify_all
#' 
#' @param annotations list, created by lc_vis
#' 
#' @export
classify_all = function(annotations) {

  # Initialize list to store all classifications
  ann_list <- list()

  for (smp in names(annotations)) {
    
    # Check whether there are actually annotations
    if (length(annotations[[smp]])>0) {

      # Initialize tibble to store classifications
      tmp_tibble <- tibble::tibble(
        barcode = annotations[[smp]][[1]]$barcode,
        classification = "none")
      
      # Get cell types that were assigned in this sample
      types_sample <- names(annotations[[smp]])
      
      for (type in types_sample) {
        
        tmp_tibble <- tmp_tibble %>%
          dplyr::mutate(tmp = annotations[[smp]][[type]]$all) %>%
          dplyr::mutate(classification = dplyr::case_when(
            # ambigious if classification is not nonce and another class is positive for the cell
            !(classification=="none") & .data$tmp ~ "ambigious",
            # if classification is nonce and the cell is positive for the class assign the type
            classification=="none" & .data$tmp ~ type,
            TRUE ~ classification
          )) }
      
      class_tmp <- tmp_tibble %>% dplyr::pull(.data$classification)
      
        ann_list[[smp]] <- tibble::tibble(
          id = annotations[[smp]][[1]]$barcode,
          classification = class_tmp
        )
      }
    }
  return(ann_list)
}

#' Function to get all rules in a given sample
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' 
#' @name get_all_rules
#' 
#' @param rules list, containing the created rules for each sample for each cell type
#' @param sample_select string, selected sample which must be in sample_info$sample
#' @param all_types logical, TRUE: return rules for all cell types in sample, FALSE: returns rules for current cell type
#' @param current_class string, currently selected cell type
#' 
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

#' Populate the app environment
#' @name initialize_app_env
#' 
#' @param app_env R environment, app enviroment
#' @param cells data.frame, containing cell metadata
#' @param counts sparse matrix, containing counts (cols = cells, genes = rows)
#' @param pc_space matrix, containing PCA coordinates (cols = coordiantes, rows = cells)
#' @param embedding matrix, containing coordinates for two-dimensional embedding (cols = coordiantes, rows = cells)
#' @param nn list, containing nearest neighbor information (created by run_nn function)
#' @param k numeric, number of neighbors considered for smoothing (default k = 50)
#' 
initialize_app_env = function(app_env, cells, counts, pc_space, embedding, nn, k) {
  app_env$cells <- cells
  app_env$counts <- counts
  app_env$pc_space <- pc_space
  app_env$embedding <- embedding
  app_env$nn <- nn
}

#' Checking whether a valid cell type was supplied by the user
#' @name check_cell_type
#' 
#' @param ct string, cell type input by user
#' @param app_env R environment, app enviroment
#' 
check_cell_type = function(ct, app_env) {
  if (!(is.character(ct))) {
    
    app_env$display <- "Enter a string"
    return(FALSE)
    
  } else if (ct == "") {
    
    app_env$display <- "Name cannot be empty"
    return(FALSE)
    
  } else {
    return(TRUE)
  }
}


