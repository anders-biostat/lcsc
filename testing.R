# Some testing
tmp$sample_info

sample_select = "PTC1"
selection <- cells$sample==sample_select

correction <- tmp$sample_info %>% dplyr::filter(sample==sample_select) %>% dplyr::pull(correction)

nn_subset <- list(idx = (nn$idx[selection,] - correction),
                          dists = nn$dist[selection,])

marker <- "Pdgfra"
gamma <- 0.5
x <- smoothed_fraction(counts, marker, selection, nn_subset, gamma)

