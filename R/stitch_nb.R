#' Compile Stitched Neighborhood Graph
#'
#' @description The main orchestrating function that sequentially chains wormhole
#' injection, absorption sweeps, and island stitching to ensure a single
#' connected spatial topology component.
#'
#' @param nb An 'nb' neighbor list object from the \code{spdep} package.
#' @param coord A 2-column matrix of projected grid centroids (e.g., UTM coordinates).
#' @param wormhole A data frame containing custom wormhole coordinates, or NULL.
#' @param max_pass Maximum absorption sweep loop parameters. Defaults to 10.
#'
#' @return A custom object of class \code{stitch_nb} containing the fully connected spatial graph.
#' @export
stitch_nb <- function(nb, coord, wormhole = NULL, max_pass = 10) {
  coord <- validate_coord(coord)

  # Stage 1: Add Wormholes
  if (!is.null(wormhole)) {
    for (i in seq_len(nrow(wormhole))) {
      nb <- add_wormhole(nb, wormhole[i, ], coord)
    }
  }

  # Stage 2: Sweep Absorption
  nb <- sweep_absorp(nb, coord, max_pass = max_pass)

  # Stage 3: Stitch Islands
  nb <- stitch_island(nb, coord)

  # Integrity Convergence Audit
  g <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE) |>
    igraph::graph_from_adjacency_matrix(mode = "undirected")
  final_comp <- igraph::components(g)$no

  if (final_comp > 1) {
    warning(paste("Graph convergence incomplete. Remaining disconnected sub-components:", final_comp))
  }

  class(nb) <- c("stitch_nb", class(nb))
  return(nb)
}

#' @param object An object of class \code{stitch_nb}.
#' @param ... Additional arguments passed to downstream methods.
#'
#' @rdname stitch_nb
#' @method summary stitch_nb
#' @export
summary.stitch_nb <- function(object, ...) {
  cat("=== Connected Graph Target Matrix (stitch_nb) ===\n")
  NextMethod("summary")
}
