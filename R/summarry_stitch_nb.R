#' Summary of a Stitched Neighborhood Graph
#' @param object An object of class `stitch_nb`.
#' @param ... Additional arguments affecting the summary produced.
#' @export
summary.stitch_nb <- function(object, ...) {
  cat("--- Stitched Neighborhood Graph ---\n")
  
  # Calculate component statistics
  g <- spdep::nb2mat(object, style = "B", zero.policy = TRUE) |>
    igraph::graph_from_adjacency_matrix(mode = "undirected")
  comp <- igraph::components(g)
  
  cat(sprintf("Total Nodes: %d\n", length(object)))
  cat(sprintf("Connected Components: %d\n", comp$no))
  cat(sprintf("Size of Giant Component: %d\n", max(comp$csize)))
  
  if (comp$no > 1) {
    cat("Warning: Graph is still disconnected!\n")
  } else {
    cat("Success: Graph is fully connected.\n")
  }
}