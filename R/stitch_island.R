#' Stitch Isolated Islands
#' @description Explicitly isolates regions with zero neighbors and forces a topological connection to their nearest mainland counterpart.
#' @param nb An 'nb' neighbor list object.
#' @param coord A 2-column matrix of projected grid coordinates (UTM).
#' @return An updated 'nb' object.
#' @export
stitch_island <- function(nb, coord) {
  coord <- validate_coord(coord)
  cards <- spdep::card(nb)
  islands <- which(cards == 0)

  if (length(islands) == 0) return(nb)

  mainland <- which(cards > 0)
  if (length(mainland) == 0) stop("No valid connected components found to snap islands to.")

  for (isl in islands) {
    nn <- RANN::nn2(coord[mainland, , drop = FALSE], coord[isl, , drop = FALSE], k = 1)
    nearest_mainland_node <- mainland[nn$nn.idx[1]]
    nb <- nb_set_edge(nb, isl, nearest_mainland_node)
  }
  return(nb)
}
