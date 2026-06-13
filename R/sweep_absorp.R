#' Sweep Absorption Pass
#' @description Iteratively evaluates disjoint sub-graph fragments, snapping them to the main giant network component.
#' @param nb An nb neighbor list object.
#' @param coord A 2-column matrix of projected grid coordinates (UTM).
#' @param max_pass Integer. Max allowed graph evaluation loops.
#' @param search_radius Numeric. Max spatial distance allowed for an absorption link.
#' @return Updated 'nb' object.
#' @export
sweep_absorp <- function(nb, coord, max_pass = 10, search_radius = 1000) {
  
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
      
      # Respect the user's search radius before snapping
      if (nn$nn.dists[best_idx] <= search_radius) {
         frag_node <- frag_sample[best_idx]
         main_node <- main_indices[nn$nn.idx[best_idx]]
         nb <- nb_set_edge(nb, frag_node, main_node)
         edges_added <- edges_added + 1
      }
    }
    
    # If a pass occurs but no fragments were close enough to absorb, break early
    if (edges_added == 0) break 
  }
  return(nb)
}