#' Iteratively replace NA values in a numeric matrix
#'
#' This function iteratively replaces missing values in a numeric matrix with
#' the mean of neighbouring cell values. Both orthogonal and diagonal neighbours
#' are included when calculating the mean value. The function is used by
#' \code{get_prop_reflected} to remove missing values from ground elevation
#' matrices.
#'
#' @param M A numeric matrix possibly containing missing values
#'
#' @return A numeric matrix with missing values replaced. If all values were
#'   missing in the input matrix, a warning message is issues and the matrix is
#'   returned unchanged.
#'
#' @export
#'
nibble_matrix <- function(M) {
  stopifnot(is.matrix(M))
  di <- -1:1

  clamp <- function(j) max(1,min(j)):min(ncol(M),max(j))

  na.inds <- which(is.na(M), arr.ind = TRUE)
  if (nrow(na.inds) == nrow(M)) {
    warning("All values are missing in input matrix")
    return(M)
  }

  while (nrow(na.inds) > 0) {
    nfixed <- 0
    for (ina in 1:nrow(na.inds)) {
      ir <- na.inds[ina,1]
      ic <- na.inds[ina,2]
      kr <- clamp(ir + di)
      kc <- clamp(ic + di)
      k.inds <- as.matrix(expand.grid(kr, kc))
      val <- mean(M[k.inds], na.rm = TRUE)
      if (!is.na(val)) {
        M[ir,ic] <- val
        nfixed <- nfixed + 1
      }
    }

    # If we failed to fix any of the NA cells on this iteration
    # it's time to give up
    if (nfixed == 0) break
    else na.inds <- which(is.na(M), arr.ind = TRUE)
  }

  M
}


#' Create a raster brick from \code{get_prop_reflected} results
#'
#' This function takes a result list object returned by
#' \code{\link{get_prop_reflected}} and creates a \code{RasterBrick} object
#' where the cell values in each layer represent the proportion of rays returned
#' for each voxel column.
#'
#' @param res A \code{TLSResult} object as returned by \code{get_prop_reflected}.
#'
#' @param elem The result element to convert to a \code{RasterBrick}. One of:
#'   'preflected' (default), 'nreflected', 'nthrough', 'nblocked', 'lower', 'upper'.
#'   An abbreviated name will work if unique. The 'lower' and 'upper' options only
#'   apply if the result object contains arrays for lower and upper bounds of the
#'   proportion of reflected rays.
#'
#' @export
#'
res_to_brick <- function(res,
                         elem = c('preflected', 'nreflected', 'nthrough', 'nblocked',
                                  'lower', 'upper')) {

  if (!inherits(res, "TLSResult")) {
    stop("Expected a list of class 'TLSResult' as returned by get_prop_reflected")
  }

  elem = match.arg(elem)

  if (elem == 'lower' || elem == 'upper') {
    ii <- grepl("^bound_", names(res))
    if (sum(ii) != 2) {
      stop("This result object does not have arrays for lower and upper bounds")
    }

    bounds.names <- names(res)[ii]

    # Assume bounds objects are listed in order lower, upper
    if (elem == 'lower') {
      ar <- res[[bounds.names[1]]]
    } else {
      ar <- res[[bounds.names[2]]]
    }
  } else { # something other than bounds
    ar <- res[[elem]]
  }

  b <- raster::brick(ar,
                     xmn = res$xylims[1], xmx = res$xylims[2],
                     ymn = res$xylims[1], ymx = res$xylims[2])

  z <- seq(res$zlims[1], res$zlims[2], by = res$zsize)[-1]
  names(b) <- paste0("to_", ifelse(z < 0, "minus_", ""), abs(z), "m")

  b
}

