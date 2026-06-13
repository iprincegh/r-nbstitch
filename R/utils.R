#' Validate and Extract Coordinate Projection Matrix
#' @description Internal helper to handle sf geometries automatically. If unprojected
#' data is provided alongside a target CRS, nbstitch does the transformation internally.
#' @param coord An 'sf' or 'sfc' object.
#' @param crs Optional. Target coordinate reference system (e.g., 3857, 32632).
#' @importFrom sf st_is_longlat st_transform st_coordinates
#' @return A 2-column numeric matrix of validated projected coordinates.
#' @keywords internal
validate_coord <- function(coord, crs = NULL) {
  if (!inherits(coord, "sf") && !inherits(coord, "sfc")) {
    stop("Input 'coord' must be an 'sf' or 'sfc' object.")
  }
  
  # Handle unprojected data automatically if a CRS is provided
  if (sf::st_is_longlat(coord)) {
    if (is.null(crs)) {
      stop("Unprojected coordinates (Lat/Lon) detected. Please provide a 'crs' argument (e.g., crs = 32632) to let nbstitch project it automatically.")
    }
    message("Detected unprojected data. Transforming to specified CRS internally...")
    coord <- sf::st_transform(coord, crs = crs)
  }
  
  # Extract coordinates as a 2-column matrix for distance calculations
  coords_mat <- sf::st_coordinates(coord)
  if (ncol(coords_mat) < 2) {
    stop("Failed to extract 2D coordinates from the sf object.")
  }
  
  return(coords_mat[, 1:2])
}

#' Symmetrically Set a Graph Edge
#' @description Adds undirected links to an nb list object.
#' @param nb An 'nb' neighbor list object.
#' @param i Integer index of first node.
#' @param j Integer index of second node.
#' @return Updated 'nb' object.
#' @keywords internal
nb_set_edge <- function(nb, i, j) {
  if (!(j %in% nb[[i]])) nb[[i]] <- sort(c(nb[[i]], j))
  if (!(i %in% nb[[j]])) nb[[j]] <- sort(c(nb[[j]], i))
  # Remove 0s if they were previously islands
  nb[[i]] <- nb[[i]][nb[[i]] > 0]
  nb[[j]] <- nb[[j]][nb[[j]] > 0]
  return(nb)
}