#include <RcppArmadillo.h>
//((Rcpp::depends(RcppArmadillo)))

using namespace Rcpp;

// Count total and (optionally) reflected rays for each voxel.
// If only counting total rays per voxel, set countReflect = FALSE and ground = NULL;
// otherwise set countReflected = TRUE and ground = matrix of ground elevations
// relative to scanner height.
//
// [[Rcpp::export]]
List cpp_count_voxel_rays(
    NumericMatrix scandata,
    int numcols, int numheights,
    double xysize, double minxy, double maxxy,
    double zsize, double minz, double maxz,
    double maxdist,
    bool countReflected = false,
    Rcpp::Nullable<NumericMatrix> groundarg = R_NilValue) {

  arma::ucube reflectedrays;
  arma::ucube throughrays;
  arma::ucube blockedrays;
  NumericMatrix ground;

  // tolerance to use for checking that x,y,z coordinates are
  // within bounds
  const double TOL = 1e-4;

  if (countReflected) {
    if (groundarg.isNull()) {
      throw std::runtime_error("Matrix of ground elevations must be provided to count reflections");
    }

    ground = NumericMatrix(groundarg);

    // check ground elevation matrix
    if (!(ground.ncol() == numcols && ground.nrow() == numcols)) {
      stop("Incorrect ground matrix size. Expected %i rows and columns", numcols);
    }

    reflectedrays.set_size(numcols, numcols, numheights);
    reflectedrays.fill(0);

    throughrays.set_size(numcols, numcols, numheights);
    throughrays.fill(0);

    blockedrays.set_size(numcols, numcols, numheights);
    blockedrays.fill(0);
  }

  arma::ucube totalrays(numcols, numcols, numheights);
  totalrays.fill(0);

  int maxtotalsegments = (int) std::trunc(maxdist / zsize) + 1;

  for (int rayindex = 0; rayindex < scandata.nrow(); rayindex++) {
    double xo = scandata(rayindex, 0);
    double yo = scandata(rayindex, 1);
    double zo = scandata(rayindex, 2);
    double dreflect = sqrt(xo*xo + yo*yo + zo*zo);

    // the points between segments can be written as
    //    x,y,z = k*delta*(xo,yo,zo)
    // where k is an integer and delta is a constant that ensures
    // each segment has length zsize
    double delta = zsize / dreflect;

    int segstodo = maxtotalsegments;

    int maxXsegments = (int) std::round(maxxy / std::abs(delta * xo));
    segstodo = std::min(segstodo, maxXsegments);

    int maxYsegments = (int) std::round(maxxy / std::abs(delta * yo));
    segstodo = std::min(segstodo, maxYsegments);

    // max Z segments
    int maxZsegments;
    if (zo < 0)
      maxZsegments = std::round(std::abs(minz / (delta * zo)));
    else
      maxZsegments = std::round(std::abs(maxz / (delta * zo)));

    segstodo = std::min(segstodo, maxZsegments);

    // loop over segments for this ray
    int reflectedSegment = std::trunc(dreflect / zsize);
    for (int k = 0; k < segstodo; k++) {
      // creating points at midpoints (k=0.5,1.5,etc)
      double xn = (k + 0.5) * delta * xo;
      double yn = (k + 0.5) * delta * yo;
      double zn = (k + 0.5) * delta * zo;

      // guard against points right on boundary
      if (xn > minxy + TOL &&
          xn < maxxy - TOL &&
          yn > minxy + TOL &&
          yn < maxxy - TOL) {

        // x col index
        int xcol = std::trunc((xn - minxy) / xysize);

        // y row index, calculated as per raster package row indexing
        // with the first row being max Y ordinate and last row being
        // min Y ordinate
        int yrow = numcols - 1 - std::trunc((yn - minxy) / xysize);

        int zlayer = std::trunc((zn - minz) / zsize);

        if (zlayer < 0 || zlayer >= numheights) {
          stop("out of bounds: zlayer (%i) for z (%f)", zlayer, zn);
        }

        // increment the total number of observed rays in the voxel
        totalrays(yrow, xcol, zlayer)++ ;

        if (countReflected) {
          double znadj = zn - ground(yrow, xcol);
          if (znadj > 0) {
            int zlayeradj = std::trunc(znadj / zsize);
            if (k < reflectedSegment) throughrays(yrow, xcol, zlayeradj)++ ;
            else if (k == reflectedSegment) reflectedrays(yrow, xcol, zlayeradj)++ ;
            else blockedrays(yrow, xcol, zlayeradj)++ ;
          }
        }
      }
    }
  }

  List res;

  // Return as a vector with elements in the right order to
  // fill a 2D matrix (heights as cols)
  res["total"] = wrap(totalrays.cols(0, numcols-1));

  if (countReflected) {
    res["through"] = wrap(throughrays.cols(0, numcols-1));
    res["reflected"] = wrap(reflectedrays.cols(0, numcols-1));
    res["blocked"] = wrap(blockedrays.cols(0, numcols-1));
  }

  return res;
}

