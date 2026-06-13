#' Compile Stitched Neighborhood Graph
#'
#' @description The main orchestrating function that sequentially chains wormhole
#' injection, absorption sweeps, and island stitching. Handles internal sf extraction
#' and projection.
#'
#' @param nb An 'nb' neighbor list object from the `spdep` package.
#' @param coord An 'sf' object containing spatial geometries.
#' @param wormhole A data frame containing custom wormhole coordinates, or NULL.
#' @param crs Integer/Character. Target CRS (e.g., EPSG code) to automatically project Lat/Lon data.
#' @param max_dist Numeric. Maximum distance allowed for an infrastructure edge. Defaults to 5000.
#' @param search_radius Numeric. Distance limit for absorption sweeps. Defaults to 1000.
#' @param k_neighbors Integer. Number of nearest neighbors to connect for isolated islands. Defaults to 1.
#' @param max_pass Integer. Maximum absorption sweep loops. Defaults to 10.
#' @param ... Additional arguments passed to internal helper functions.
#'
#' @return A custom object of class `stitch_nb` containing the fully connected spatial graph.
#' @export
stitch_nb <- function(nb, 
                      coord, 
                      wormhole = NULL, 
                      crs = NULL,
                      max_dist = 5000,
                      search_radius = 1000,
                      k_neighbors = 1,
                      max_pass = 10,
                      ...) {
  
  # The gatekeeper: Let nbstitch handle the sf transformation and extraction
  coord_mat <- validate_coord(coord, crs = crs)

  # Stage 1: Add Wormholes
  if (!is.null(wormhole)) {
    for (i in seq_len(nrow(wormhole))) {
      nb <- add_wormhole(nb, wormhole[i, ], coord_mat, max_dist = max_dist)
    }
  }

  # Stage 2: Sweep Absorption
  nb <- sweep_absorp(nb, coord_mat, max_pass = max_pass, search_radius = search_radius)

  # Stage 3: Stitch Islands
  nb <- stitch_island(nb, coord_mat, k_neighbors = k_neighbors)

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