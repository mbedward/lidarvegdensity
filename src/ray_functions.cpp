#include <RcppArmadillo.h>
//((Rcpp::depends(RcppArmadillo)))

using namespace Rcpp;

// Count number of rays passing through each voxel
//
// [[Rcpp::export]]
NumericVector rays_in_voxels(NumericMatrix scandata,
                             int numcols, int numheights,
                             double xysize, double minxy, double maxxy,
                             double zsize, double minz, double maxz,
                             double maxdist) {

  arma::ucube prismrays(numcols, numcols, numheights);
  prismrays.fill(0);

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

    // loop over segments
    for (int k = 0; k < segstodo; k++) {
      // creating points at midpoints (k=0.5,1.5,etc)
      double xn = (k + 0.5) * delta * xo;
      double yn = (k + 0.5) * delta * yo;
      double zn = (k + 0.5) * delta * zo;

      // need to know which column the midpoint is in
      int xcol = std::trunc((xn - minxy) / xysize);
      int ycol = std::trunc((yn - minxy) / xysize);

      // and need the z columns
      int zcol = std::trunc((zn - minz) / zsize);

      // increment the total number of rays in the prism
      prismrays(xcol, ycol, zcol)++ ;
    }
  }

  // Return as a vector with elements in the right order to
  // fill a 2D matrix (heights as cols)
  return wrap(prismrays.cols(0, numcols-1));
}

