#include <R.h>
#include <Rinternals.h>
#define NCOL 6

// functions that parse paf files. To start with a function for parsing the difference
// string only. Ideally this should be done as part of parsing the whole file, but
// I don't want to write the whole parser at the moment.


// ONLY HANDLES SHORT FORM OF THE STRING!!!
// and does not yet handle introns; but we should be able to sort
// that out afterwards with a bit more complicated stuff
SEXP parse_diff_string(SEXP diff_str_r, SEXP rstart_r, SEXP qstart_r, SEXP qend_r, 
		       SEXP fwd_strand_r, SEXP include_strings_r)
{
  if(TYPEOF(diff_str_r) != STRSXP || length(diff_str_r) < 1)
    error("The argument should be a non-null character vector");
  int n = length(diff_str_r);
  if(TYPEOF(rstart_r) != INTSXP || TYPEOF(qstart_r) != INTSXP || TYPEOF(qend_r) != INTSXP)
    error("rstart, qstart and qend should be integer vectors");
  if(TYPEOF(fwd_strand_r) != LGLSXP || TYPEOF(include_strings_r) != LGLSXP)
    error("fwd_strand and include_strings should be logical vectors");
  if(length(rstart_r) != n || length(qstart_r) != n || length(qend_r) != n || length(fwd_strand_r) != n)
    error("cs, rstart and rend need to be the same length");
  // we will need a list
  int *rstart = INTEGER(rstart_r);
  int *qstart = INTEGER(qstart_r);
  int *qend = INTEGER(qend_r);
  int *fwd_strand = LOGICAL(fwd_strand_r);
  int include_strings = LOGICAL(include_strings_r)[0];
  SEXP ret_data = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT( ret_data, 0, allocVector(VECSXP, n));
  SEXP alignments = VECTOR_ELT(ret_data, 0);
  if(include_strings)
    SET_VECTOR_ELT( ret_data, 1, allocVector(VECSXP, n));
  SEXP mut_strings_list = VECTOR_ELT(ret_data, 1);
  //  SEXP ret_data = PROTECT(allocVector(VECSXP, n));
  // Each entry will be converted to a matrix with the following columnames
  //  const int ncol = 6;
  const char *colnames[NCOL] = {"op", "length", "mut", "ref.pos", "query.o", "query.pos"};
  SEXP dim_names = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dim_names, 0, NILSXP);
  SET_VECTOR_ELT(dim_names, 1, allocVector(STRSXP, NCOL));
  /* SEXP colnames_r = PROTECT(allocVector(STRSXP, 5)); */
  for(int i=0; i < NCOL; ++i)
    SET_STRING_ELT( VECTOR_ELT(dim_names, 1), i, mkChar(colnames[i]));

  for(int i=0; i < n; ++i){
    int op_n=0;
    const char *str = CHAR(STRING_ELT(diff_str_r, i));
    const char *p = str;
    while(*p){
      if( *p == ':' || *p == '*' || *p == '+' || *p == '-' )
	++op_n;
      ++p;
    }
    // allocate useful memory. Unfortunately all counts need to be in ints; and to use
    // in R we probably want to do something clever to convert letters to numbers.
    // We can do a series of bitwise operations to convert
    // : -> 0  * -> 2  + -> 3  -  -> 1
    // and the same bitwise operation gives unique numbers if we include the last
    // three bits.. 
    // but use the numeric value of the operator instead: it is much simpler..
    // because we can then use utf8ToInt and intToUtf8 to convert to strings..
    // return a matrix with columns:
    // 1. operation (numeric value)
    // 2. length associated with the operation
    // 3. mutation associated encoded in the second and first 8 bits of the integer
    // 4. reference position (always offset + start)
    // 5. query offset
    // 6. query position
    // increment op_n by 1 to have space for the final row, which should always
    // be a ':' operation (since we do not support the the long format).
    op_n++;
    SET_VECTOR_ELT( alignments, i, allocMatrix(INTSXP, op_n, NCOL) );
    SET_VECTOR_ELT( dim_names, 0, allocVector(STRSXP, op_n));
    setAttrib( VECTOR_ELT(alignments, i), R_DimNamesSymbol, dim_names );

    if(include_strings)
      SET_VECTOR_ELT( mut_strings_list, i, allocVector(STRSXP, op_n) );
    SEXP mut_strings = include_strings ? VECTOR_ELT(mut_strings_list, i) : NILSXP;
    
    int *op_matrix = INTEGER( VECTOR_ELT( alignments, i ));
    int *ops = op_matrix;
    int *length = op_matrix + op_n;
    int *mut = op_matrix + op_n * 2;
    int mut_default = (int)'.' | ( (int)'.' << 8);
    int *ref = op_matrix + op_n * 3;
    int *query_o = op_matrix + op_n * 4;
    int *query_p = op_matrix + op_n * 5;
    memset( (void*)op_matrix, 0, op_n * NCOL * sizeof(int));
    int row = 0;
    p = str;
    char int_buf[19];  // for holding numbers use atoi to convert to integer
    // set these so they count from 1 based offset
    int ref_offset = rstart[i] + 1;
    int query_offset = fwd_strand[i] ? 1 : 0; // fwd_strand[i] ? qstart[i] : qend[i];
    int l = 0;
    while(*p && row < (op_n-1)){
       mut[row] = mut_default;
       if(*p == ':'){  // aligned and matching
	++p;
	int dig = 0;
	while(*p >= '0' && *p <= '9' && dig < sizeof(int_buf)-1){
	  int_buf[dig] = *p;
	  ++dig;
	  ++p;
	}
	int_buf[dig] = 0;
	l = atoi(int_buf); // overflow error?
	int_buf[dig] = 0;
	ops[row] = (int)':';
	length[row] = l;
	ref[row] = ref_offset;
	query_o[row] = query_offset;
	query_p[row] = fwd_strand[i] ? qstart[i] + query_offset : qend[i] - query_offset;
	ref_offset += l;
	query_offset += l;
	++row;
	continue;
      }
      if(*p == '*'){  // SNP
	unsigned int r_c = (unsigned int)p[1]; // could be 0
	unsigned int q_c = (p[1] != 0) ? (unsigned int)p[2] : 0;
	if(include_strings && q_c)
	  SET_STRING_ELT( mut_strings, row, mkCharLen(p+1, 2) );
	ops[row] = (int)'*';
	length[row] = 1;
	mut[row] = (r_c << 8) | q_c;
	ref[row] = ref_offset;
	query_o[row] = query_offset;
	query_p[row] = fwd_strand[i] ? qstart[i] + query_offset : qend[i] - query_offset;
	ref_offset++;
	query_offset++;
	p += (r_c) ? 2 : 1;
	p += (q_c) ? 1 : 0;
	++row;
	continue;
      }
      if(*p == '+'){  // Insertion into the query
	int l = 0;
	++p;
	const char *beg = p;
	while(*p && (*p >= 'a' && *p <= 'z')){
	  ++p;
	  ++l;
	}
	if(include_strings)
	  SET_STRING_ELT(mut_strings, row, mkCharLen(beg, l));
	ops[row] = (int)'+';
	length[row] = l;
	ref[row] = ref_offset;
	query_o[row] = query_offset;
	query_p[row] = fwd_strand[i] ? qstart[i] + query_offset : qend[i] - query_offset;
	query_offset += l;
	++row;
	continue;
      }
      if(*p == '-'){ // Deletion from query
	int l = 0;
	++p;
	const char *beg = p;
	while(*p && (*p >= 'a' && *p <= 'z')){
	  ++p;
	  ++l;
	}
	if(include_strings)
	  SET_STRING_ELT(mut_strings, row, mkCharLen(beg, l));
	ops[row] = (int)'-';
	length[row] = l;
	ref[row] = ref_offset;
	query_o[row] = query_offset;
	query_p[row] = fwd_strand[i] ? qstart[i] + query_offset : qend[i] - query_offset;
	ref_offset += l;
	++row;
	continue;
      }
      // since we don't support = and ~ at the moment, that should be everything we need to do
      // I may also want to extract the indels into sequences, but first I'm more interested in their
      // lengths.
    }
    // set the last row (i.e. the end of the last operation)
    ops[row] = (int)':';
    length[row] = l;
    mut[row] = mut_default;
    ref[row] = ref_offset;
    query_o[row] = query_offset;
    query_p[row] = fwd_strand[i] ? qstart[i] + query_offset : qend[i] - query_offset;
  }
  UNPROTECT(2);
  return(ret_data);
}
