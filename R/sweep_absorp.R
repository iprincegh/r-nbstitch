#' Sweep Absorption Pass
#' @description Iteratively evaluates disjoint sub-graph fragments, snapping them to the main giant network component.
#' @param nb An nb neighbor list object.
#' @param coord A 2-column matrix of projected grid coordinates (UTM).
#' @param max_pass Max allowed graph evaluation loops. Defaults to 10.
#' @return Updated 'nb' object.
#' @export
sweep_absorp <- function(nb, coord, max_pass = 10) {
  coord <- validate_coord(coord)

  for (pass in seq_len(max_pass)) {
    g <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE) |>
      igraph::graph_from_adjacency_matrix(mode = "undirected")
    comp <- igraph::components(g)

    if (comp$no <= 1) break

    comp_sizes <- comp$csize
    giant_id <- which.max(comp_sizes)
    main_indices <- which(comp$membership == giant_id)

    edges_added <- 0

    for (cid in seq_len(comp$no)) {
      if (cid == giant_id) next

      frag_indices <- which(comp$membership == cid)
      sample_size <- min(5, length(frag_indices))
      frag_sample <- frag_indices[seq(1, length(frag_indices), length.out = sample_size)]

      nn <- RANN::nn2(coord[main_indices, , drop = FALSE], coord[frag_sample, , drop = FALSE], k = 1)
      best_idx <- which.min(nn$nn.dists)

      node_frag <- frag_sample[best_idx]
      node_giant <- main_indices[nn$nn.idx[best_idx]]

      if (!(node_giant %in% nb[[node_frag]])) {
        nb <- nb_set_edge(nb, node_frag, node_giant)
        edges_added <- edges_added + 1
      }
    }
    if (edges_added == 0) break
  }
  return(nb)
}
