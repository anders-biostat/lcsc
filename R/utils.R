#' Making life a little easier
#' @importFrom dplyr %>%
#' @export
dplyr::`%>%`

#' Generate nearest neighbor per sample
#' @name run_nn
#' @importFrom dplyr %>%
#' @param cells data.frame with cell metadata
#' @param pc_space PCA coordinates
#' @param k number of neighbors considered for smoothing
#' @param dim number of PC dimensions used for nearest neighbor graph
#' @export
run_nn <- function(cells, pc_space, k=50, dim=30){

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
#' @export
run_classification = function(annotations) {

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
