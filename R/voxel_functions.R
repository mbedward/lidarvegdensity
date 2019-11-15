#' Estimate the potential maximum number of returns for voxels
#'
#' This function takes data from a control scan that was performed in a closed
#' room (i.e. no gaps or reflective surfaces) where all laser pulses are
#' expected to provide a return.
#'
#' @param scandata A three column matrix of XYZ coordinates for return positions
#'   relative to the scanner.
#'
#' @param xysize The horizontal size (i.e. width in the X-Y plane) of voxels in
#'   metres.
#'
#' @param zsize The vertical size (height) of voxels in metres.
#'
#' @param numcols The integer number of voxel colums in the X and Y direction.
#'   Ashcroft et al. recommend that this be an odd value so that the scanner
#'   lies in the centre column, but this is not enforced here.
#'
#' @param zmax The maximum height of the scan (relative to the scanner) in
#'   metres. Default value is 35.25 as used by Ashcroft et al. This will usually
#'   be a positive value (uppermost voxel above scanner) but zero and negative
#'   values are also valid.
#'
#' @param zmin The minimum height of the scan (relative to the scanner) in
#'   metres. Must be less than \code{zmax}. Defaults to \code{-zmax} if zmax is
#'   positive. If \code{zmax} is zero or negative an explicit value must be
#'   provided.
#'
#' @return A matrix with voxel horizontal position as rows (X varies fastest,
#'   then Y) and voxel height as columns.
#'
#' @export
#'
estimate_max_returns <- function(scandata, xysize, zsize, numcols,
                                 zmax = 35.25, zmin = NULL) {

  xysize <- .gtzero( .single_numeric(xysize) )
  zsize <- .gtzero( .single_numeric(zsize) )
  numcols <- .gtzero( .single_numeric(numcols) )

  zmax <- .single_numeric(zmax)

  if (is.null(zmin)) {
    if (zmax > 0) zmin <- -zmax
    else stop("When zmax <= 0 you must provide a value for zmin")
  } else {
    zmin <- .single_numeric(zmin)
    if (zmin >= zmax) stop("Value of zmin must be less than zmax")
  }

  # total limits of the columns
  xylims <- c(-1, 1) * numcols * xysize / 2

  # max and min height in metres are defaulted to 35.25m and -35.25m
  # ensure they are multiples of zsize
  maxz <- (trunc(zmax / zsize) + 1) * zsize
  minz <- (trunc(zmin / zsize) - 1) * zsize
  numheights <- (maxz - minz) / zsize

  # calculate maximum distance of points from scanner
  maxzdist <- max(abs(c(maxz, minz)))
  maxdist <- sqrt(2*xylims[2]^2 + maxzdist^2)

  counts <- rays_in_voxels(
    scandata,
    numcols, numheights,
    xysize, xylims[1], xylims[2],
    zsize, minz, maxz,
    maxdist
  )

  matrix(counts, ncol = numheights)
}


# R reference implementation (no Rcpp) for testing
#
R_estimate_max_returns <- function(scandata, xysize, zsize, numcols,
                                   zmax = 35.25, zmin = NULL) {

  xysize <- .gtzero( .single_numeric(xysize) )
  zsize <- .gtzero( .single_numeric(zsize) )
  numcols <- .gtzero( .single_numeric(numcols) )

  zmax <- .single_numeric(zmax)

  if (is.null(zmin)) {
    if (zmax > 0) zmin <- -zmax
    else stop("When zmax <= 0 you must provide a value for zmin")
  } else {
    zmin <- .single_numeric(zmin)
    if (zmin >= zmax) stop("Value of zmin must be less than zmax")
  }

  if (numcols %% 2 == 0) {
    numcols <- numcols + 1
    message("Adjusted numcols to ", numcols)
  }

  # total limits of the columns
  xylims <- c(-1, 1) * numcols * xysize / 2

  # max and min height in metres are defaulted to 35.25m and -35.25m
  # ensure they are multiples of zsize
  maxz <- (trunc(zmax / zsize) + 1) * zsize
  minz <- (trunc(zmin / zsize) - 1) * zsize
  numheights <- (maxz - minz) / zsize

  # calculate maximum distance of points from scanner
  maxzdist <- max(abs(c(maxz, minz)))
  maxdist <- sqrt(2*xylims[2]^2 + maxzdist^2)

  # create array to store number of ray segments in each rectangular prism
  prismrays <- array(0, dim = c(numcols, numcols, numheights))

  # We are going to go through point by point, creating a laser ray and dividing it
  # into segments that are each zsize in length. This is equivalent to
  # the vertical distance between each height class.

  # iterate through each ray
  pb <- txtProgressBar(max = nrow(scandata), style = 3)
  for (rayindex in 1:nrow(scandata)) {
    if (rayindex %% 1e4 == 0) setTxtProgressBar(pb, rayindex)

    # get the ray data for this point
    xo <- scandata[rayindex,1]
    yo <- scandata[rayindex,2]
    zo <- scandata[rayindex,3]
    dreflect <- sqrt(xo^2 + yo^2 + zo^2)

    # the points between segments can be written as x,y,z = k*delta*(xo,yo,zo) where k is
    # an integer and delta is a constant that ensures each segment has length zsize
    delta <- zsize/dreflect

    # work out how many segments within range
    maxtotsegments <- trunc(maxdist / zsize) + 1
    maxxsegments <- round(xylims[2] / abs(delta * xo))
    maxysegments <- round(xylims[2] / abs(delta * yo))

    if (zo < 0)
      maxzsegments <- round(abs(minz / (delta * zo)))
    else
      maxzsegments <- round(abs(maxz / (delta * zo)))

    numsegments <- min(maxtotsegments, maxxsegments, maxysegments, maxzsegments)

    # loop over segments
    for (k in 1:numsegments) {
      # creating points at midpoints (k=0.5,1.5,etc)
      xn = (k - 0.5) * delta * xo
      yn = (k - 0.5) * delta * yo
      zn = (k - 0.5) * delta * zo

      # need to know which column the midpoint is in
      xcol <- trunc((xn - xylims[1]) / xysize) + 1
      ycol <- trunc((yn - xylims[1]) / xysize) + 1

      # and need the z columns
      zcol <- trunc((zn - minz) / zsize) + 1

      # inc total number of rays in the prism
      prismrays[xcol, ycol, zcol] <- prismrays[xcol, ycol, zcol] + 1
    }
  }
  setTxtProgressBar(pb, nrow(scandata))
  close(pb)

  browser()

  # Flatten 3D array to matrix with heights as cols
  matrix(c(prismrays), nrow = numcols * numcols, ncol = numheights)
}


.single_numeric <- function(x) {
  if (is.numeric(x) && length(x) == 1) {
    x
  } else {
    nm <- deparse(substitute(x))
    paste(nm, "should be a single numeric value")
  }
}

.single_integer <- function(x) {
  if (is.integer(x) && length(x) == 1) {
    x
  } else {
    nm <- deparse(substitute(x))
    paste(nm, "should be a single integer value")
  }
}

.gtzero <- function(x) {
  if (is.numeric(x) && all(x > 0)) {
    x
  } else {
    nm <- deparse(substitute(x))
    paste(nm, "should be greater than zero")
  }
}


