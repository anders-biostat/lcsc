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

x <- get_all_rules(rules, "PTC1", TRUE, tmp$current_class)

g_env <- globalenv()
app_env <- tmp

g_env$annotations[[app_env$sample_select]][[app_env$current_class]]

# Resetting the `all` column
g_env$annotations[[app_env$sample_select]][[app_env$current_class]]$all <- rep(NULL, length(app_env$classification))

g_env$annotations[[app_env$sample_select]][[app_env$current_class]][app_env$marker] <- app_env$classification

c_names <- colnames(g_env$annotations[[app_env$sample_select]][[app_env$current_class]])


system.time(
# For checking whether a cell adheres to all rules we have to check all column (=genes) are TRUE
x <- g_env$annotations[[app_env$sample_select]][[app_env$current_class]] %>%
  dplyr::rowwise() %>%
  dplyr::mutate(all = all(dplyr::c_across(cols = setdiff(c_names, c("barcode", "all"))))) # select all columns but barcode & all
)


system.time(
  tmp <- g_env$annotations[[app_env$sample_select]][[app_env$current_class]] %>%
    select(setdiff(c_names, c("barcode", "all")))
)

system.time(
  p <- apply(tmp, MARGIN=1, FUN=all)
)


p <- t( g_env$annotations[[app_env$sample_select]][[app_env$current_class]] )

y <- g_env$annotations[[app_env$sample_select]][[app_env$current_class]]

Pdgrf = sample(c(TRUE, FALSE), size=nrow(g_env$annotations[[app_env$sample_select]][[app_env$current_class]]), replace=TRUE)
tmp <- g_env$annotations[[app_env$sample_select]][[app_env$current_class]] %>%
  select(setdiff(c_names, c("barcode", "all"))) %>%
  cbind(Pdgrf)

p <- apply(tmp, MARGIN=1, FUN=all)
