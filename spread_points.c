#include <R.h>
#include <Rinternals.h>
#include <string.h>

// This takes a single vector of values indicating positions along
// an axis and provides positions along a perpendicular axis which
// will spread out the points in a reasonable manner.
// Note that the values must be sorted before use

// This function should be included in general drawing functions.

double dist(double x1, double y1, double x2, double y2){
  return( sqrt((x1-x2) * (x1-x2) + (y1-y2) * (y1-y2)) );
}


SEXP spread_points( SEXP y_r, SEXP min_dist_r ){
  if(!isVector(y_r) || !isVector(min_dist_r))
    error("Both arguments should be vectors");
  if(!isReal(y_r) || !isReal(min_dist_r))
    error("Both arguments should be real values");
  if(!length(y_r) || !length(min_dist_r))
    error("Both arguments must contain values");
  // y should be sorted
  if(isUnsorted( y_r, FALSE ))
    error("y values must be sorted");

  int yn = length(y_r);
  double *y = REAL(y_r);
  double min_dist = asReal( min_dist_r );

  // we will return a set of offset positions.
  // refer to the offsets as the x-position for convenience
  SEXP x_r = PROTECT(allocVector( REALSXP, yn ));
  double *x = REAL( x_r ); 
  memset( x, 0, sizeof(double) * yn );
  
  // We always accept positions that are below
  for(int i=1; i < yn; ++i){
    // check every previous one to see if we need to do soemthing
    double x_shift = min_dist / 3;
    double x1 = 0, x2 = 0;
    double d1, d2;
    d1 = d2 = 0;
    while( (d1 < min_dist && d2 < min_dist) ){
      d1 = d2 = 2 * min_dist;
      int j=i;
      while( j > 0 ){
	--j;
	if( y[i] - y[j] > min_dist )
	  break;
	d1 = d1 < min_dist ? d1 : dist(x1, y[i], x[j], y[j]);
	d2 = d2 < min_dist ? d2 : dist(x2, y[i], x[j], y[j]);
	if( d1 < min_dist && d2 < min_dist )
	  break;
      }
      if( d1 >= min_dist && fabs(x1) <= fabs(x2) )
	x[i] = x1;
      if( d2 >= min_dist && fabs(x2) <= fabs(x1) )
	x[i] = x2;
      x1 = d1 < min_dist ? x1 + x_shift : x1;
      x2 = d2 < min_dist ? x2 - x_shift : x2;
    }
  }
  UNPROTECT(1);
  return( x_r );
}
