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
#' @return A named list with the following elements:
#'   \describe{
#'   \item{returns}{A matrix of maximum return values for voxels. Rows
#'   correspond to horizontal position, ordered (ascending) by X position within
#'   Y position. Columns correspond to vertical position, ordered (ascending) by
#'   height.}
#'   \item{xysize}{The horizontal voxel width.}
#'   \item{xlims}{A vector of min and max values for X and Y ordinates.}
#'   \item{zsize}{The vertical voxel width.}
#'   \item{zlims}{A vector or min and max values for Z (height).}
#'   \item{numcols}{NUmber of voxels in each horizontal dimension.}
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
  zlims <- c('min' = (trunc(zmin / zsize) - 1) * zsize,
             'max' = (trunc(zmax / zsize) + 1) * zsize)

  numheights <- diff(zlims) / zsize
  stopifnot(numheights > 0)

  # calculate maximum distance of points from scanner
  maxzdist <- max(abs(c(maxz, minz)))
  maxdist <- sqrt(2*xylims[2]^2 + maxzdist^2)

  # Call Rcpp function to count returns for each voxel
  counts <- rays_in_voxels(
    scandata,
    numcols, numheights,
    xysize, xylims[1], xylims[2],
    zsize, zlims[1], zlims[2],
    maxdist
  )

  # Format as matrix where rows are horizontal voxel position ordered
  # by X within Y; columns are vertical position ordered by Z
  counts <- matrix(counts, ncol = numheights)

  list(returns = counts,
       xysize = xysize, xylims = xylims,
       zsize = zsize, zlims = zlims,
       numcols = numcols, numheights = numheights)
}


#' Main function to process scans
#'
#' This function takes scan data for a vegetation site, together with a matrix of
#' max returns per voxel generated using the function \code{get_max_returns()},
#' and calculates something...
#'
#' @param scandata A three column matrix of XYZ coordinates for return positions
#'   relative to the scanner.
#'
#' @param refdata A named list, as returned by the functcion
#'   \code{get_max_returns()}, containing a matrix of the estimated maximum
#'   number of returns for each voxel plus the voxel and scan dimensions.
#'
#' @param ground A square matrix where values are ground elevation for each
#'   voxel. The dimensions of the matrix must match the value of element
#'   \code{numcols} in \code{refdata}.
#'
#' @return A matrix where cell values are proportion of rays reflected for each voxel
#'   based on the reference data for the potential maximum number.
#'
#'
# rawscandata is the scan you are interested in, nullmodel is the array produced NULL12,
# ground is a DTM at the same resolution as columnsize/heightsize. I have included
# comments on what is going on below. My comments are not aligned with the rest of
# the text

get_prop_reflected <- function(scandata, refdata, ground) {

  xysize <- refdata$xysize
  zsize <- refdata$zsize

  if (nrow(ground) != refdata$numcols ||
      ncol(ground) != refdata$numcols) {

    stop("The matrix of ground elevations should be square, with number of\n",
         "rows and columns equal to the 'numcols' value in refdata")
  }

  # Determine the height range and number of voxel height levels.
  # We clamp the height limits to multiples of zsize.
  maxz <- (trunc(max(scandata[,3]) / zsize) + 1) * zsize
  minz <- (trunc(min(ground, na.rm = TRUE) / zsize) - 1) * zsize
  numheights <- (maxz - minz) / zsize

  # Height offset value used to relate vertical voxel position in scan to
  # position in refdata
  # TODO - is this really necessary?
  zoffset <- (minz + 150) / zsize

  # maximum distance based on the processing extent defined earlier.
  # calculate maximum distance of points from scanner
  maxzdist <- max(abs(maxz),abs(minz))
  maxdist <- sqrt(2*maxxy^2 + maxzdist^2)

  #These are the arrays where the data is stored
  # create arrays to store number of ray segments in each rectangular prism

  # prismrays stores according to raw z values (minz to maxz) and is the total
  # number of rays (blocked before reaching this prism, reflected from this prism +
  # passed through this prism)

  # throughrays and returnrays are stored according to ground level (0 to maxz-minz)
  # throughrays is number that pass through the prism without reflecting back
  # reflectedrays is number that are reflected back from within this prism
  numprisms <- numcolumns * numcolumns * numheights
  prismrays <- rep(0, numprisms)
  dim(prismrays) <- c(numcolumns, numcolumns, numheights)
  throughrays <- reflectedrays <- prismrays
  print(paste("Finished creating output arrays: ",date()))

  # We are going to go through point by point, creating a laser ray and dividing it
  # into segments that are each zsize in length. This is equivalent to
  # the vertical distance between each height class.
  # For each segment we will determine if the ray is blocked (reflection occurred before
  # segment), reflected (reflection occurred in segment) or through (reflection
  # occurred after segment) and add the totals to our output arrays


  #This is the main process for the function and takes the longets. It runs through every beam and finds which voxels it passed through and which voxel it recived a return from. Given there are over 14 million returns it can take a long time to process. There are reasonably good details about what is going on here.
  # iterate through each ray
  for (rayindex in 1:length(scandata[,1]))
  {
    debug <- FALSE
    # provide debug info every 50000 points
    # comment out these lines to remove debug
    #if ((rayindex-1) == (trunc((rayindex-1)/50000)*50000))
    #{
    #	debug <- TRUE
    #}

    # get the ray data for this point
    xo <- scandata[rayindex,1]
    yo <- scandata[rayindex,2]
    zo <- scandata[rayindex,3]
    lengthtoreflection <- sqrt(xo^2 + yo^2 + zo^2)

    # the points between segments can be written as x,y,z = k.c.(xo,yo,zo) where k is
    # an integer and c is a constant that ensures each segment has length zsize
    c <- zsize/lengthtoreflection

    # work out how many segments within range
    maxtotsegments <- trunc(maxdist / zsize) + 1
    maxxsegments <- trunc(maxxy / abs(c * xo) + 0.4999999999)
    maxysegments <- trunc(maxxy / abs(c * yo) + 0.4999999999)
    if (zo < 0)
    {
      maxzsegments <- trunc(abs(minz / (c * zo)) + 0.4999999999)
    }		else
    {
      maxzsegments <- trunc(abs(maxz / (c * zo)) + 0.4999999999)
    }
    numsegments <- min(maxtotsegments, maxxsegments, maxysegments, maxzsegments)

    # work out the segment the reflection is in
    reflectedsegment <- trunc(lengthtoreflection / zsize) + 1

    if (debug)
    {
      print(paste("raynum: ",rayindex, " ", date()))
      print(paste("ray: ",xo, " ", yo, " ", zo))
      print(paste("lengths: ",lengthtoreflection, " ", c))
      print(paste("maxs: ",maxtotsegments," ", maxxsegments, " ",maxysegments," ",
                  maxzsegments))
      print(paste("reflected segment: ",reflectedsegment))
    }

    # loop over segments
    for (k in 1:numsegments)
    {

      # creating points at midpoints (k=0.5,1.5,etc)
      xn = (k - 0.5) * c * xo
      yn = (k - 0.5) * c * yo
      zn = (k - 0.5) * c * zo

      # need to know which column the midpoint is in
      xcol <- trunc((xn - minxy) / columnsize) + 1
      ycol <- trunc((yn - minxy) / columnsize) + 1

      # and need the z columns (raw for prismrays, adj for others)
      zcol <- trunc((zn - minz) / zsize) + 1

      if (debug)
      {
        print(paste(k, " ",xn," ",yn," ",zn," ",xcol," ",ycol," ",zcol))
        print(ground[xcol,ycol])
      }
      # and need the adjusted height / col
      znadj <- zn - ground[xcol,ycol]
      zcoladj <- trunc(znadj / zsize)  + 1

      if (debug)
      {
        print(paste(znadj," ",zcoladj))
      }

      # inc total number of rays in the prism
      prismrays[xcol, ycol, zcol] <- prismrays[xcol, ycol, zcol] + 1

      # and add to blocked etc, if it is above ground level
      if (znadj > 0)
      {
        if (k == reflectedsegment)
        {
          if (debug)
          {
            print("Reflected")
          }
          reflectedrays[xcol,ycol,zcoladj] <- reflectedrays[xcol,ycol,zcoladj] + 1
        }
        else if (k < reflectedsegment)
        {
          if (debug)
          {
            print("Through")
          }
          throughrays[xcol,ycol,zcoladj] <- throughrays[xcol,ycol,zcoladj] + 1
        }
      }
      # we won't count blocked rays, as canopy cover is determined by ratio of reflected
      # to reflected + through. Blocked didn't get far enough to provide any indication
      # of canopy cover, but we still need to loop through all segments to increment
      # prismrays
    }
  }

  print(paste("Finished processing laser points: ",date()))

  #This section also takes a while, but not as long.

  # We've gone through all the points provided in the scan, but some rays went right through the canopy without
  # being returned. We need to add these rays into the 'through' tally by comaprison with a null model in
  # a closed room.
  for (xindex in 1:numcolumns)
  {
    for (yindex in 1:numcolumns)
    {
      for (zindex in 1:numheights)
      {
        zminlimit <- (zindex-1) * zsize + minz
        zmaxlimit <- zindex * zsize + minz
        if (prismrays[xindex,yindex,zindex] < refdata[xindex+(yindex-1)*numcolumns,zindex+zoffset])
        {
          # we have less rays in this prism than expected. Let's add them in
          # at random.
          for (index in 1:(refdata[xindex+(yindex-1)*numcolumns,zindex+zoffset]
                           - prismrays[xindex,yindex,zindex]))
          {
            # Create a random point in this prism
            # for x and y we can use original indexes, but z will change
            # according to ground level calculation
            #z <- zminlimit + runif(1) * (zmaxlimit - zminlimit)
            z <- zminlimit + 0.5 * (zmaxlimit - zminlimit)
            zadj <- z - ground[xindex,yindex]
            zcoladj <- trunc(zadj / zsize)  + 1
            throughrays[xindex,yindex,zcoladj] <-
              throughrays[xindex,yindex,zcoladj] + 1
          }
        }
      }
    }
  }

  print(paste("Finished adding random through lasers based on null model: ",date()))

  # create an array to hold our output. Number of rows = numcolumns x numcolumns.
  # number of columns = numheights. Output is ratio of reflected:(reflect + through)
  output<- rep(0, (numcolumns * numcolumns * numheights))
  dim(output) <- c(numcolumns * numcolumns, numheights)

  # Fill in the data
  for (xindex in 1:numcolumns)
  {
    for (yindex in 1:numcolumns)
    {
      for (zindex in 1:numheights)
      {
        output[xindex+(yindex-1)*numcolumns,zindex] <-
          reflectedrays[xindex,yindex,zindex] /
          (reflectedrays[xindex,yindex,zindex] +
             throughrays[xindex,yindex,zindex])
      }
    }
  }
  print(paste("Finished processing: ",date()))

  return(output)
} # endfunction



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




