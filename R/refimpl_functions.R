# R reference implementation (no Rcpp) for testing
#
R_get_max_returns <- function(scandata, xysize, zsize, numcols,
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

  #if (numcols %% 2 == 0) {
  #  numcols <- numcols + 1
  #  message("Adjusted numcols to ", numcols)
  #}

  # total limits of the columns
  xylims <- c(-1, 1) * numcols * xysize / 2

  # max and min height in metres are defaulted to 35.25m and -35.25m
  # ensure they are multiples of zsize
  zlims.working <- c(
    (trunc(zmin / zsize) - 1) * zsize,
    (trunc(zmax / zsize) + 1) * zsize
  )
  numheights <- abs(diff(zlims.working) / zsize)

  # calculate maximum distance of points from scanner
  maxzdist <- max(abs(zlims.working))
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
      maxzsegments <- round(abs(zlims.working[1] / (delta * zo)))
    else
      maxzsegments <- round(abs(zlims.working[2] / (delta * zo)))

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
      zcol <- trunc((zn - zlims.working[1]) / zsize) + 1

      # inc total number of rays in the prism
      prismrays[xcol, ycol, zcol] <- prismrays[xcol, ycol, zcol] + 1
    }
  }
  setTxtProgressBar(pb, nrow(scandata))
  close(pb)

  list(returns = prismrays,
       xysize = xysize, xylims = xylims,
       zsize = zsize, zlims = c(zmin, zmax),
       numcols = numcols, numheights = numheights)
}



# R reference implementation (no Rcpp) for testing
#
# ground - either a square matrix with the same number of rows and columns as
# the refdata object; or a single value. If a matrix, values are ground elevations
# relative to scanner. A single value is equivalent to a matrix of uniform values.
#
# Returns a list with meta-data elements and a 3D array where the first two
# dimensions (rows and columns) match the input refdata object and the third
# (vertical) dimension represents voxel layers from the minimum ground height
# to the maximum scan height.
#
R_get_prop_reflected <- function(scandata, refdata, ground) {

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

  #These are the arrays where the data is stored
  # create arrays to store number of ray segments in each rectangular prism

  # prismrays stores according to raw z values (minz to maxz) and is the total
  # number of rays (blocked before reaching this prism, reflected from this prism +
  # passed through this prism)

  # throughrays and returnrays are stored according to ground level (0 to maxz-minz)
  # throughrays is number that pass through the prism without reflecting back
  # reflectedrays is number that are reflected back from within this prism
  prismrays <- array(0, dim = c(numcols, numcols, numheights))
  throughrays <- reflectedrays <- prismrays

  # We are going to go through point by point, creating a laser ray and dividing it
  # into segments that are each zsize in length. This is equivalent to
  # the vertical distance between each height class.
  # For each segment we will determine if the ray is blocked (reflection occurred before
  # segment), reflected (reflection occurred in segment) or through (reflection
  # occurred after segment) and add the totals to our output arrays


  # This is the main process for the function and takes the longest. It runs through
  # every beam and finds which voxels it passed through and which voxel it recived
  # a return from. Given there are over 14 million returns it can take a long time
  # to process. There are reasonably good details about what is going on here.
  # iterate through each ray
  for (rayindex in 1:nrow(scandata)) {
    xo <- scandata[rayindex,1]
    yo <- scandata[rayindex,2]
    zo <- scandata[rayindex,3]
    dreflect <- sqrt(xo^2 + yo^2 + zo^2)

    # the points between segments can be written as x,y,z = k*delta*(xo,yo,zo) where k is
    # an integer and delta is a constant that ensures each segment has length zsize
    delta <- zsize / dreflect

    # work out how many segments within range
    maxtotsegments <- trunc(maxdist / zsize) + 1
    maxxsegments <- round(xylims[2] / abs(delta * xo))
    maxysegments <- round(xylims[2] / abs(delta * yo))

    if (zo < 0)
      maxzsegments <- round(abs(minz / (delta * zo)))
    else
      maxzsegments <- round(abs(maxz / (delta * zo)))

    numsegments <- min(maxtotsegments, maxxsegments, maxysegments, maxzsegments)

    # work out the segment the reflection is in
    reflectedsegment <- trunc(dreflect / zsize) + 1


    # loop over segments
    for (k in 1:numsegments) {
      # creating points at midpoints (k=0.5,1.5,etc)
      xn = (k - 0.5) * delta * xo
      yn = (k - 0.5) * delta * yo
      zn = (k - 0.5) * delta * zo

      # need to know which column the midpoint is in
      xcol <- trunc((xn - xylims[1]) / xysize) + 1
      ycol <- trunc((yn - xylims[1]) / xysize) + 1

      # and need the z columns (raw for prismrays, adj for others)
      zcol <- trunc((zn - minz) / zsize) + 1

      znadj <- zn - ground[xcol, ycol]
      zcoladj <- trunc(znadj / zsize) + 1

      # inc total number of rays in the prism
      prismrays[xcol, ycol, zcol] <- prismrays[xcol, ycol, zcol] + 1

      # and add to reflected and through if it is above ground level
      if (znadj > 0) {
        if (k == reflectedsegment)
          reflectedrays[xcol,ycol,zcoladj] <- reflectedrays[xcol,ycol,zcoladj] + 1
        else if (k < reflectedsegment)
          throughrays[xcol,ycol,zcoladj] <- throughrays[xcol,ycol,zcoladj] + 1
        #else
        #  blockedrays[xcol,ycol,zcol] <- blockedrays[xcol,ycol,zcol] + 1
      }
    }
  }

  # We've gone through all the points provided in the scan, but some rays went
  # right through the vegetation without being returned. We need to add these rays
  # into the 'through' tally by comaprison with a null model in a closed room.
  for (xindex in 1:numcols) {
    for (yindex in 1:numcols) {
      for (zindex in 1:numheights) {
        # Check if we have less rays in the voxel than expected (because
        # some went through). If so, increase the through count.
        expected <- refdata$returns[xindex, yindex, zindex + ref.zoffset]
        if (prismrays[xindex,yindex,zindex] < expected) {
          nadd <- expected - prismrays[xindex, yindex, zindex]

          zminlimit <- (zindex-1) * zsize + minz
          zmaxlimit <- zindex * zsize + minz

          z <- zminlimit + 0.5 * (zmaxlimit - zminlimit)
          zadj <- z - ground[xindex,yindex]
          zcoladj <- trunc(zadj / zsize)  + 1
          throughrays[xindex,yindex,zcoladj] <- throughrays[xindex,yindex,zcoladj] + nadd
        }
      }
    }
  }

  # Return proportions
  preflect <- reflectedrays / (reflectedrays + throughrays)

  # Guard against NaN which occur when a voxel in both reflectedrays
  # and throughrays has a zero count
  preflect[is.nan(preflect)] <- 0

  list(preflect = preflect,
       xysize = xysize, xylims = xylims,
       zsize = zsize, zlims = c(minz, maxz),
       numcols = numcols, numheights = numheights,
       totalrays = prismrays,
       reflectedrays = reflectedrays,
       throughrays = throughrays)
}


