#' Stitch Nodes via Correlated Temporal Activity
#' @description Links isolated sub-graphs to the main spatial component by comparing
#' temporal activity profiles (e.g., mobile-phone visit volumes).
#' @param nb An 'nb' neighbor list object.
#' @param activity A numeric matrix of temporal activity (rows = spatial units, columns = time).
#' @param cor_threshold Numeric. Minimum Pearson correlation required to add an edge. Defaults to 0.6.
#' @param n_candidates Integer. Number of most active cells in the main component to search against. Defaults to 300.
#' @return An updated 'nb' object.
#' @export
stitch_activity <- function(nb, activity, cor_threshold = 0.6, n_candidates = 300) {
  if (nrow(activity) != length(nb)) {
    stop("'activity' matrix rows must exactly match the length of 'nb'.")
  }

  g <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE) |>
    igraph::graph_from_adjacency_matrix(mode = "undirected")
  comps <- igraph::components(g)

  if (comps$no == 1) return(nb)

  main_id <- which.max(comps$csize)
  main_nodes <- which(comps$membership == main_id)
  isolate_nodes <- which(comps$membership != main_id)

  main_activity_sums <- rowSums(activity[main_nodes, , drop = FALSE], na.rm = TRUE)
  top_main_nodes <- main_nodes[order(main_activity_sums, decreasing = TRUE)[1:min(n_candidates, length(main_nodes))]]

  edges_added <- 0

  for (i in isolate_nodes) {
    profile_i <- activity[i, ]
    if (all(is.na(profile_i) | profile_i == 0)) next

    cors <- apply(activity[top_main_nodes, , drop = FALSE], 1, function(x) {
      if (stats::sd(x, na.rm = TRUE) == 0 || stats::sd(profile_i, na.rm = TRUE) == 0) return(-1)
      stats::cor(profile_i, x, use = "pairwise.complete.obs")
    })

    max_cor <- max(cors, na.rm = TRUE)

    if (max_cor >= cor_threshold) {
      best_match <- top_main_nodes[which.max(cors)]
      nb <- nb_set_edge(nb, i, best_match)
      edges_added <- edges_added + 1
    }
  }

  message(sprintf("Activity Stitching: Added %d edges (r >= %s).", edges_added, cor_threshold))
  return(nb)
}