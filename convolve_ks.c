#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>  // for abs()
#include <string.h> // for memset()

// convolve a vector of values with a symmetric kernel
// The user should specify a vector of values of length n 
// and half of the kernel of length m
// m should be smaller than n and should not be even
// We consider the first point of the kernel as the center and
// define it's radius (r) as m-1
//
// The index of the kernel (0..(m-1)) specifies the distance
// from the values it is multiplied such that the return
// value 
// v.conv[i] = sum(j=i-r:i+r ; kern[k] * v[j] ) / sum(kern[abs(i-j)])
//
// where k = abs( j - i )
// i.e. the distance from the kernel center.
// note that sum(kern) refers to the part of the kernel used and 
// will in most kases be: 2 * kern[1:r] + kern[0]
// This implies a complexity of m x n
// which is trivial for short kernels.

// All values will be assumed to be doubles for convenience

// pass the return data in by reference so that its memory can be allocated by
// R functions
// val_conv should be zeroed; but that is really up to the user
// Note that we are using int as an offset. This is because that is what length( )
// returns. It also makes the abs( ) call and subtractions more easy to deal with
// But it is kind of dirty.
void convolve_values( double *val, int v_n, double *kern, int k_n, double *val_conv )
{
  int r = k_n - 1;
  for(int i=0; i < v_n; ++i){
    int beg = i - r > 0 ? i - r : 0;
    int end = i + r < v_n ? i + r : v_n - 1; // inclusive range
    double k_sum = 0;
    //    Rprintf("%d :  %f   %f\n", i, val[i], val_conv[i]);
    for(int j=beg; j <= end; ++j){
      int ko = abs(i-j);
      k_sum += kern[ ko ];
      val_conv[i] += val[j] * kern[ko];
      //      Rprintf("\t|%d %d %f  %f|", j, ko, val[j], val_conv[i]);
    }
    val_conv[i] /= k_sum;
    //    Rprintf("\n\t->%f\n", val_conv[i]);
  }
}

// ks is short for kernel, symmetric
// both val_r and kern_r should be positive length
// double vectors
SEXP convolve_ks(SEXP val_r, SEXP kern_r)
{
  if(!isVector(val_r) || !isVector(kern_r) || !isReal(val_r) || !isReal(kern_r))
    error("Both arguments should be vectors of doubles");
  
  int v_n = length(val_r);
  int k_n = length(kern_r);
  if(k_n < 1 || k_n % 2 == 0)
    error("The kernel should have a positive odd length greater than 0");
  if(v_n < k_n)
    error("The kernel must be smaller than the value vector");
  
  // allocate the return vector
  SEXP convolved_r = PROTECT( allocVector(REALSXP, v_n) );
  double *convolved = REAL(convolved_r);
  memset( convolved, 0, sizeof(double) * v_n );

  convolve_values( REAL(val_r), v_n, REAL(kern_r), k_n, convolved );
  UNPROTECT(1);
  return( convolved_r );
}

