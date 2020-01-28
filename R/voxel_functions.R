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
#'   metres. Must be less than \code{zmax}. Defaults to minus \code{zmax} if
#'   zmax is positive. If \code{zmax} is zero or negative an explicit value must
#'   be provided.
#'
#' @return A named list with the following elements:
#'   \describe{
#'   \item{returns}{An array (dimensions: X, Y, Z) of maximum expected number
#'     of rays passing through each voxel.}
#'   \item{xysize}{The horizontal voxel width.}
#'   \item{xlims}{A vector of min and max values for X and Y ordinates.}
#'   \item{zsize}{The vertical voxel width.}
#'   \item{zlims}{A vector or min and max values for Z (height).}
#'   \item{numcols}{Number of voxels in each horizontal dimension.}
#'   \item{numheights}{Number of voxels in vertical dimension.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Estimate returns for regular cubic voxels with 0.25m side length,
#' # for a scan block with horizontal dimensions 60x60 voxels
#' refdata <- get_max_returns(control.scan,
#'                            xysize = 0.25, zsize = 0.25, numcols = 60,
#'                            zmax = 100)
#' }
#'
#' @export
#'
get_max_returns <- function(scandata, xysize, zsize, numcols,
                            zmax = 35.25, zmin = NULL) {

  xysize <- .gtzero( .single_numeric(xysize) )
  zsize <- .gtzero( .single_numeric(zsize) )
  numcols <- .gtzero( .single_numeric(numcols) )

  xylims <- c(-1, 1) * numcols * xysize / 2
  names(xylims) <- c('min', 'max')

  zmax <- .single_numeric(zmax)

  if (is.null(zmin)) {
    if (zmax > 0) zmin <- -zmax
    else stop("When zmax <= 0 you must provide a value for zmin")
  } else {
    zmin <- .single_numeric(zmin)
    if (zmin >= zmax) stop("Value of zmin must be less than zmax")
  }

  # max and min height in metres, ensuring that they are
  # multiples of zsize
  zlims.working <- c(
    (trunc(zmin / zsize) - 1) * zsize,
    (trunc(zmax / zsize) + 1) * zsize)

  numheights <- abs(diff(zlims.working) / zsize)
  stopifnot(numheights > 0)

  # calculate maximum distance of points from scanner
  maxzdist <- max(abs(zlims.working))
  maxdist <- sqrt(2*xylims[2]^2 + maxzdist^2)

  # Call Rcpp function to count the number of rays passing through
  # each voxel
  counts <- cpp_count_voxel_rays(
    scandata,
    numcols, numheights,
    xysize, xylims[1], xylims[2],
    zsize, zlims.working[1], zlims.working[2],
    maxdist
  )

  # Reshape the returned vector into an array
  # (the cpp_count_voxel_rays function should have ordered the
  # vector appropriately)
  counts <- array(counts$total, dim = c(numcols, numcols, numheights))

  list(returns = counts,
       xysize = xysize, xylims = xylims,
       zsize = zsize, zlims = c(zmin, zmax),
       numcols = numcols, numheights = numheights)
}


#' Main function to process scans
#'
#' This function takes scan data for a vegetation site, together with a reference
#' object giving the maximum expected number of returns per voxel, and calculates
#' the proportion of rays reflected for each voxel. The reference object is
#' generated using the function \code{get_max_returns()}. Optionally, Bayesian
#' credible intervals can be calculated for the proportion values by specifying
#' one or more probability values via the \code{probs} argument. See details below.
#'
#' Bayesian credible intervals on proportions are based on a non-informative
#' Jeffrey's prior distribution. For a given probability p, this involves taking
#' the lower (\code{(1-p)/2}) and upper (\code{(1+p)/2}) quantiles from a
#' \code{Beta(1/2 + r, 1/2 + t)} distribution, where r is the number of
#' reflected rays (returns) for the voxel, and t is the number of rays passing
#' through the voxel without being reflected. For voxels where no rays were
#' reflected, the lower bound is set to 0. For voxels where all rays were
#' reflected, the upper bound is set to 1. If a voxel was not crossed by any
#' rays, its calculated proportion will be missing (\code{NA}) and its interval
#' will also be set to \code{NA}.
#'
#' @param scandata Either a three column matrix or data frame of XYZ coordinates
#'   for return positions relative to the scanner; or a \code{LAS} object as
#'   returned by \code{lidR::readLAS}.
#'
#' @param refdata A named list, as returned by the function
#'   \code{get_max_returns()}, containing an array of the estimated maximum
#'   number of returns for each voxel plus the voxel and scan dimensions.
#'
#' @param ground A square matrix where values are ground elevation for each
#'   voxel. The dimensions of the matrix must match the value of element
#'   \code{numcols} in \code{refdata}.
#'
#' @param probs Optionally, one or more probability values (0 to 1). If
#'   provided, bounds on the proportion of reflected rays are calculated for
#'   each voxel. For example, with \code{probs = c(0.5, 0.9)} lower and upper
#'   50 and 90 percent bounds will be calculated.
#'
#' @return A named list with the following elements:
#'   \describe{
#'   \item{preflect}{An array (dimensions: X, Y, Z) of the proportion of rays
#'     reflected for each voxel.}
#'   \item{bounds_XX}{If argument probs was specified, arrays for lower and
#'     upper bounds on calculated proportions will be included,
#'     e.g. setting probs = 0.9 will results in two arrays 'bounds_0.05'
#'     and 'bounds_0.95'}
#'   \item{xysize}{The horizontal voxel width.}
#'   \item{xlims}{A vector of min and max values for X and Y ordinates.}
#'   \item{zsize}{The vertical voxel width.}
#'   \item{zlims}{A vector or min and max values for Z (height).}
#'   \item{numcols}{NUmber of voxels in each horizontal dimension.}
#'   \item{numheights}{Number of voxels in vertical dimension.}
#'   }
#'
#' @export
#'
get_prop_reflected <- function(scandata, refdata, ground, probs = NULL) {

  if (inherits(scandata, "LAS")) {
    scandata <- as.matrix(scandata@data[, c("X", "Y", "Z")])

  } else if (is.matrix(scandata) || is.data.frame(scandata)) {
    if (!ncol(scandata) == 3) {
      stop("Argument scandata should be a three column matrix or data frame")
    }
    scandata <- as.matrix(scandata)
  }

  # Check probs argument
  if (length(probs) > 0) {
    if (!is.numeric(probs) || !is.vector(probs) || !all(probs > 0 & probs < 1)) {
      stop("Argument probs should be a vector of one or more values between 0 and 1")
    }
  }

  numcols <- refdata$numcols
  xylims <- refdata$xylims
  xysize <- refdata$xysize
  zsize <- refdata$zsize

  # Check that scan points are within the bounds of the reference model
  inside.ref <- scandata[,1] >= xylims[1] &&
    scandata[,1] <= xylims[2] &&
    scandata[,2] >= xylims[1] &&
    scandata[,2] <= xylims[2]

  if (any(!inside.ref)) {
    n <- sum(!inside.ref)
    msg <- ifelse(n == 1, "point is outside", "points are outside")
    msg <- glue::glue("{n} {msg} the bounds of the reference model and will be ignored")
    warning(msg)
  }

  # Check ground elevation data
  if (is.matrix(ground)) {
    if (nrow(ground) != refdata$numcols ||
        ncol(ground) != refdata$numcols) {

      stop("The matrix of ground elevations should be square, with number of\n",
           "rows and columns equal to the 'numcols' value in refdata")
    }

  } else if (is.numeric(ground)) {
    ground <- .single_numeric(ground)
    ground <- matrix(ground, numcols, numcols)

  } else {
    stop("Argument ground should either be a square matrix with number of\n",
         "rows and columns equal to the 'numcols' value in refdata, or a\n",
         "single value for ground elevation to apply to all voxels.")
  }

  # If there are any NA values in ground, set these to the mean
  # non-NA neighbour value. If any cells have only missing neighbour values
  # set these (arbitrarily) to the min ground value and issue a warning.
  if (anyNA(ground)) {
    warn.nas <- FALSE
    di <- -1:1
    km <- matrix(0, nrow=9, ncol=2)
    clamp <- function(x) max(1,min(x)):min(numcols,max(x))
    for (ir in 1:numcols) {
      for (ic in 1:numcols) {
        x <- 0
        n <- 0
        for (kr in clamp(ir + d)) {
          for (kc in clamp(ic + d)) {
            g <- ground[kr, kc]
            if (!is.na(g)) {
              x <- x + g
              n <- n + 1
            }
          }
        }
        if (n > 0) ground[ir, ic] <- x/n
        else {
          ground[ir, ic] <- ground.minz
          warn.nas <- TRUE
        }
      }
    }

    if (warn.nas)
      warning("Ground elevation for one or more cells could not be determined.\n",
              "Cell values have been set to overall minimum value.")
  }

  # Determine the height range and number of voxel height levels.
  # We clamp the height limits to multiples of zsize.
  ground.minz <- min(ground, na.rm = TRUE)
  minz <- (trunc(ground.minz / zsize) - 1) * zsize
  maxz <- (trunc(max(scandata[,3]) / zsize) + 1) * zsize

  stopifnot(maxz > minz)
  numheights <- (maxz - minz) / zsize


  # Find the index of the height dimension in the reference data array
  # that corresponds with the ground elevation for each voxel column
  ground.minzindex <- trunc((ground - refdata$zlims[1]) / zsize) + 1

  # Z-index offset for reference model.
  ref.zoffset <- (minz - refdata$zlims[1]) / zsize

  # calculate maximum distance of points from scanner
  maxzdist <- max(abs(c(maxz, minz)))
  maxdist <- sqrt(2*xylims[2]^2 + maxzdist^2)

  # Call Rcpp function to count the total number of rays entering each voxel,
  # the number reflected, and the number passing through.
  counts <- cpp_count_voxel_rays(
    scandata,
    numcols, numheights,
    xysize, xylims[1], xylims[2],
    zsize, minz, maxz,
    maxdist,
    TRUE,
    ground
  )

  totalrays <- array(counts$total, dim = c(numcols, numcols, numheights))
  reflectedrays <- array(counts$reflected, dim = c(numcols, numcols, numheights))
  throughrays <- array(counts$through, dim = c(numcols, numcols, numheights))

  # To account for rays that went through a voxel without being returned, compare
  # right through the vegetation without being returned. We need to add these rays
  # into the 'through' tally by comaprison with a null model in a closed room.
  for (xindex in 1:numcols) {
    for (yindex in 1:numcols) {
      for (zindex in 1:numheights) {
        # Check if we have less rays in the voxel than expected (because
        # some went through). If so, increase the through count.
        expected <- refdata$returns[xindex, yindex, zindex + ref.zoffset]
        if (totalrays[xindex,yindex,zindex] < expected) {
          nadd <- expected - totalrays[xindex, yindex, zindex]

          zminlimit <- (zindex-1) * zsize + minz
          z <- zminlimit + zsize/2
          zadj <- z - ground[xindex,yindex]
          zcoladj <- trunc(zadj / zsize) + 1
          throughrays[xindex,yindex,zcoladj] <- throughrays[xindex,yindex,zcoladj] + nadd
        }
      }
    }
  }

  # Proportion of reflected rays per voxel
  preflect <- reflectedrays / (reflectedrays + throughrays)

  # Convert any NaN values, which occur when a voxel had no rays
  # to NA
  preflect[is.nan(preflect)] <- NA

  res <- list(preflect = preflect)

  # Calculate bounds if probs were specified
  if (length(probs) > 0) {
    my <- matrix(c(reflectedrays, throughrays), ncol = 2)
    qs <- sort(unique(c((1-probs)/2, (1+probs)/2)))

    bounds <- lapply(1:length(qs), function(i) {
      x <- qbeta(qs[i], 1/2 + my[,1], 1/2 + my[,2])

      if (i <= length(probs)) { # lower bound
        x[my[,1] == 0] <- 0
      } else { # upper bound
        x[my[,2] == 0] <- 1
      }

      x <- array(x, dim = dim(reflectedrays))
      x[is.na(preflect)] <- NA
      x
    })

    names(bounds) <- paste0("bound_", qs)

    res <- c(res, bounds)
  }

  # Return results
  c(res,
    list(xysize = xysize, xylims = xylims,
         zsize = zsize, zlims = c(minz, maxz),
         numcols = numcols, numheights = numheights) )
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




