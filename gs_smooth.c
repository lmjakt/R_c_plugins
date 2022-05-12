#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for memset

// Smooth a set of points in the y dimension (2)
// The points are separated in x and the smoothing
// uses a Gaussian blur. The function takes:
// 1. Positions in x
// 2. Positions in y
// 3. A pointer to the smoothened points allocated by the caller
// 4. The length of x, y, and y_v (use int for compatibility with R)
// 5. The stanadard deviation (c) of the gaussian
// 6. A minimal value beyond which points are excluded from
//    the smoothing.

// smoothes points by:
// exp( -(d^2) / (2c^2) )
// where d is the distance between points

// The length of x, y, and y_sm must be the same
// uses int as offset. I feel dirty. 
void smooth_line(double *x, double *y, double *y_sm,
	       int l, double sd, double min_f){
  for(int i=0; i < l; ++i){
    double f_sum = 0; // the sum of the smoothing factors used
    double f;
    double d;
    y_sm[i] = 0;
    for(int j=i; j > 0; --j){
      d = x[i] - x[j];
      if( (f = exp( -(d*d) / (2*sd*sd))) < min_f )
	break;
      y_sm[i] += y[j] * f;
      f_sum += f;
    }
    for(int j=i+1; j < l; ++j){
      d = x[i] - x[j];
      if( (f = exp( -(d*d) / (2*sd*sd))) < min_f )
	break;
      y_sm[i] += y[j] * f;
      f_sum += f;
    }
    y_sm[i] = y_sm[i] / f_sum;
  }
}

SEXP gs_smooth(SEXP coords_r, SEXP sd_r, SEXP min_f_r){
  if(!isMatrix(coords_r) || !isReal(coords_r))
    error("coords_r should be a real matrix");
  if(!isReal(sd_r) || !isReal(min_f_r) || length(sd_r) != 1 || length(min_f_r) != 1)
    error("Both sd_r and min_f_r should be reals of length 1");
  double sd = REAL(sd_r)[0];
  double min_f = REAL(min_f_r)[0];
  SEXP dims_r = PROTECT( getAttrib( coords_r, R_DimSymbol ));
  if(length(dims_r) != 2)
    error("coords_r should have two dimensions");
  int *dims = INTEGER(dims_r);
  if(dims[1] != 2 || dims[0] < 1)
    error("coords_r should two columns and at least one row");
  double *coords = REAL(coords_r);
  int l = dims[0];
  double *x = coords;
  double *y = coords + l;
  SEXP y_sm_r = PROTECT(allocVector(REALSXP, l));
  double *y_sm = REAL(y_sm_r);
  smooth_line( x, y, y_sm, l, sd, min_f );
  UNPROTECT(2);
  return(y_sm_r);
}
