#' Add a Manual Wormhole Edge
#' @description Maps single path anchors into logical neighborhood paths based on closest centroids.
#' @param nb An nb neighbor list object.
#' @param wormhole A 1-row data frame/tibble containing projected metrics: `from_x`, `from_y`, `to_x`, `to_y`.
#' @param coord A 2-column matrix of projected grid coordinates (UTM).
#' @return An updated 'nb' object with the added edge.
#' @export
add_wormhole <- function(nb, wormhole, coord) {
  coord <- validate_coord(coord)

  # Pinpoint closest index to starting anchor
  nn_from <- RANN::nn2(coord, matrix(c(wormhole$from_x, wormhole$from_y), ncol = 2), k = 1)
  idx_from <- nn_from$nn.idx[1]

  # Pinpoint closest index to ending anchor
  nn_to <- RANN::nn2(coord, matrix(c(wormhole$to_x, wormhole$to_y), ncol = 2), k = 1)
  idx_to <- nn_to$nn.idx[1]

  if (idx_from != idx_to) {
    nb <- nb_set_edge(nb, idx_from, idx_to)
  }
  return(nb)
}
