#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for memset
#include <string.h> // for memcpy

// segment the regions of a matrix using a (reverse)
// watershed that finds peaks in values by following
// gradients. Such peaks and regions associated with
// the peaks need to have values above some minimum.

// I am using int here rather than size_t; this is not ideal,
// but it means that I can memcpy from a vector of offsets to
// a normal R object. 
struct offset_vector {
  size_t *offsets;
  size_t length;
  size_t capacity;
};

struct offset_vector ov_init(size_t cap){
  struct offset_vector ov;
  ov.offsets = malloc( sizeof(size_t) * cap );
  ov.capacity = cap;
  ov.length = 0;
  return(ov);
}

void ov_clear(struct offset_vector *ov){
  free( ov->offsets );
  ov->offsets = 0;
  ov->capacity = 0;
  ov->length = 0;
}

// convert to 0-based rows and columns
// that can be used in R.
// note that coords should have a size of 2 * nrow * ncol
void ov_to_matrix(struct offset_vector *ov, int *coords, int nrow){
  for(size_t i=0; i < ov->length; ++i){
    coords[i] = ov->offsets[i] % nrow;
    coords[i + ov->length] = ov->offsets[i] / nrow;
  }
}

void ov_push(struct offset_vector *ov, size_t v){
  if(ov->capacity == 0){
    ov->capacity = 128;
    ov->offsets = malloc(sizeof(size_t) * ov->capacity);
    ov->length = 0;
  }
  if(ov->length == ov->capacity){
    ov->capacity = ov->capacity * 2;
    size_t *n_offsets = malloc(sizeof(size_t) * ov->capacity);
    memcpy( (void*)n_offsets, (void*)ov->offsets, sizeof(size_t) * ov->length );
    free( ov->offsets );
    ov->offsets = n_offsets;
  }
  ov->offsets[ov->length] = v;
  ov->length++;
}

// The most natural way to code this is to use recursion.
// This can result in a stack overflow though if there
// the recusion is too deep.

// recursive function that sets the identity of a given
// position in the matrix
// max_id points to the maximum id which has been assigned so far.
int assign_id(int row, int column, 
	      double *m, int *id,
	      int nrow, int ncol,
	      int* max_id, double min_v,
	      struct offset_vector *peaks){
  size_t offset = column * nrow + row;
  if(id[offset] > 0)
    return(id[offset]);
  if(m[offset] < min_v)
    return(0);
  int pos_id = 0;
  // determine the search space based on the dimensions of m
  int left = column > 0 ? column - 1 : 0;
  int right = column < (ncol-1) ? column + 1 : column;
  int bottom = row > 0 ? row - 1 : 0;
  int top = row < (nrow-1) ? row + 1 : row;
  // the positions of the maximum value..
  int max_row = -1;
  int max_column = -1;
  double max_value = m[offset];
  for(int c=left; c <= right; ++c){
    for(int r=bottom; r <= top; ++r){
      size_t o = c * nrow + r;
      if( m[o] > max_value ){
	max_row=r;
	max_column=c;
	max_value=m[o];
      }
    }
  }
  // max_row and max_column will be 0 if none of the neighbours have a higher value
  // in that case we want to increment the max_id, set the value and return the id.
  // otherwise we recurse and set the id to the return value..
  if(max_row >= 0){
    // if the identity for that has been assigned, then we take that:
    pos_id = id[ max_column * nrow + max_row ];
    if(!pos_id)
      pos_id = assign_id(max_row, max_column, m, id, nrow, ncol, max_id, min_v, peaks);
    id[offset] = pos_id;
    return(pos_id);
  }
  // if max_row is negative, then the current position is a peak
  // increment *max_id, assign and return.
  (*max_id)++;
  id[offset] = *max_id;
  ov_push( peaks, offset );
  return(*max_id);
}

void trace_ridges(int row, int column, int nrow, int ncol,
		  double *m, int *ridges, double min_v, struct offset_vector *ridge_lines){
  size_t offset = column * nrow + row;
  // if a ridge has been set (1) or masked (-1), then return
  if( ridges[offset] )
    return;
  ridges[offset] = 1;
  ov_push(ridge_lines, offset);
  // the positions of the maximum value..
  int max_i = -1;
  double max_value = min_v;
  int c[8] = {-1, -1, -1, 0, 1, 1, 1, 0};
  int r[8] = {-1, 0, 1, 1, 1, 0, -1, -1};
  for(int i=0; i < 8; ++i){
    int cc = column + c[i];
    int rr = row + r[i];
    if(cc > 0 && cc < ncol && rr > 0 && rr < nrow){
      size_t o = cc * nrow + rr;
      if( m[o] > max_value && !ridges[o] ){
	max_i = i;
	max_value=m[o];
      }
    }
  }
  if(max_i >= 0){
    //ridges[ offset + c[max_i] * nrow + r[max_i] ] = 1;
    // mask the neighbours, but do not set the ridge as we do that at the
    // beginningo of the function.
    // we could consider masking from -2:+2 rather than only -1. This would
    // reduce angles.. 
    ridges[ offset + c[(max_i+1) % 8] * nrow + r[(max_i+1) % 8] ] = -1;
    ridges[ offset + c[(max_i+2) % 8] * nrow + r[(max_i+2) % 8] ] = -1;
    ridges[ offset + c[(8+max_i-1) % 8] * nrow + r[(8+max_i-1) % 8] ] = -1;
    ridges[ offset + c[(8+max_i-2) % 8] * nrow + r[(8+max_i-2) % 8] ] = -1;
    trace_ridges( row + r[max_i], column + c[max_i], nrow, ncol,
		  m, ridges, min_v, ridge_lines );
  }
}

SEXP ws_segment(SEXP mat_r, SEXP min_v_r){
  if(!isMatrix(mat_r) || !isReal(mat_r))
    error("The first argument should a double matrix");
  if(!isVector(min_v_r) || !isReal(min_v_r) || length(min_v_r) != 1)
    error("The second argument should be a single double value");
  double min_v = REAL(min_v_r)[0];
  SEXP r_dims = PROTECT( getAttrib( mat_r, R_DimSymbol ));
  if(length(r_dims) != 2)
    error("mat_r should have two dimensions (i.e. be a matrix)");
  int *dims = INTEGER(r_dims);
  int nrow = dims[0];
  int ncol = dims[1];
  if(nrow < 1 || ncol < 1){
    UNPROTECT(1);
    error("The matrix must have non-null dimensions");
  }
  double *mat = REAL(mat_r);
  // return a list of two elements; the id matrix and the positions and peaks.
  SEXP ret_data = PROTECT( allocVector(VECSXP, 2) );
  // set up an integer matrix;
  SET_VECTOR_ELT( ret_data, 0, allocMatrix(INTSXP, nrow, ncol) );
  int *id = INTEGER( VECTOR_ELT(ret_data, 0) );
  memset((void*)id, 0, sizeof(int) * nrow * ncol );
  int max_id = 0;
  // not strictly necessary to assign here, as it should already have been done.
  struct offset_vector peaks = ov_init(32);
  for(int r=0; r < nrow; ++r){
    for(int c=0; c < ncol; ++c){
      assign_id(r, c, mat, id, nrow, ncol, &max_id, min_v, &peaks);
    }
  }
  SET_VECTOR_ELT( ret_data, 1, allocMatrix(INTSXP, peaks.length, 2) );
  ov_to_matrix( &peaks, INTEGER(VECTOR_ELT(ret_data, 1)), nrow );
  ov_clear( &peaks );
  UNPROTECT(2);
  return(ret_data);
}

SEXP ws_trace_ridges(SEXP mat_r, SEXP peaks_r, SEXP min_v_r){
  if(!isMatrix(mat_r) || !isReal(mat_r))
    error("The first argument should a double matrix");
  if(!isMatrix(peaks_r) || !isInteger(peaks_r))
    error("The second argument should be an integer matrix");
  if(!isVector(min_v_r) || !isReal(min_v_r) || length(min_v_r) != 1)
    error("The second argument should be a single double value");
  double min_v = REAL(min_v_r)[0];
  SEXP mat_dims_r = PROTECT( getAttrib( mat_r, R_DimSymbol ));
  SEXP peak_dims_r = PROTECT( getAttrib( peaks_r, R_DimSymbol ));
  if(length(mat_dims_r) != 2 || length(peak_dims_r) !=2)
    error("both mat_r and peaks_r should have two dimensions");
  int *mat_dims = INTEGER(mat_dims_r);
  int *peak_dims = INTEGER(peak_dims_r);
  int nrow=mat_dims[0];
  int ncol=mat_dims[1];
  int n_peaks = peak_dims[0];
  if(nrow < 1 || ncol < 1){
    UNPROTECT(2);
    error("The matrix must have non-null dimensions");
  }
  if(n_peaks < 1 || peak_dims[1] != 2){
    UNPROTECT(2);
    error("Peaks must be given as matrix with two columns and at least one row %d , %d", n_peaks, peak_dims[2]);
  }
  double *mat = REAL(mat_r);
  int *peaks = INTEGER(peaks_r);
  // create a trace matrix
  SEXP ret_data = PROTECT( allocVector(VECSXP, 2) );
  SET_VECTOR_ELT( ret_data, 0, allocMatrix(INTSXP, nrow, ncol));
  //  SEXP ridges_r = PROTECT(
  int *ridges = INTEGER(VECTOR_ELT(ret_data, 0));
  memset((void*)ridges, 0, sizeof(int) * nrow * ncol);
  struct offset_vector ridge_lines = ov_init(256);
  for(int i=0; i < n_peaks; ++i){
    trace_ridges(peaks[i], peaks[i+n_peaks], nrow, ncol, mat, ridges, min_v, &ridge_lines);
  }
  SET_VECTOR_ELT( ret_data, 1, allocMatrix(INTSXP, ridge_lines.length, 2) );
  ov_to_matrix( &ridge_lines, INTEGER(VECTOR_ELT(ret_data, 1)), nrow );
  ov_clear(&ridge_lines);
  UNPROTECT(3);
  return(ret_data);
}
