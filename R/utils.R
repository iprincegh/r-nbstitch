#' Validate Coordinate Projection Matrix or sf Object
#' @description Internal helper to enforce projected coordinate metrics (e.g., UTM),
#' check for geodetic unprojected reference systems, and extract 2D coordinate
#' matrices from native sf spatial structures automatically.
#' @param coord A 2-column numeric matrix or an 'sf' or 'sfc' object.
#' @return A 2-column numeric matrix of validated projected coordinates.
#' @keywords internal
validate_coord <- function(coord) {
  # If an sf object is provided, verify projection and extract coordinate matrix
  if (inherits(coord, c("sf", "sfc"))) {
    if (sf::st_is_longlat(coord)) {
      stop("Geodetic coordinates (Lat/Lon) detected in sf object. 'nbstitch' requires a projected coordinate reference system (e.g., UTM).")
    }
    coord <- sf::st_coordinates(coord)[, 1:2, drop = FALSE]
  }

  if (!is.matrix(coord) || ncol(coord) != 2) {
    stop("Argument 'coord' must be a 2-column matrix or an 'sf' object containing projected coordinates.")
  }

  if (max(abs(coord[, 1])) <= 180 && max(abs(coord[, 2])) <= 90) {
    stop("Geodetic coordinates (Lat/Lon) detected. 'nbstitch' requires projected coordinates (e.g., UTM) where Euclidean distances are valid.")
  }

  return(coord)
}

#' Symmetrically Set a Graph Edge
#' @description Adds undirected links to an nb list object.
#' @param nb An 'nb' neighbor list object.
#' @param i Integer index of first node.
#' @param j Integer index of second node.
#' @keywords internal
nb_set_edge <- function(nb, i, j) {
  nb[[i]] <- sort(unique(c(nb[[i]], as.integer(j))))
  nb[[j]] <- sort(unique(c(nb[[j]], as.integer(i))))
  return(nb)
}
