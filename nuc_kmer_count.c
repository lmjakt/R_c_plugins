#include <R.h>
#include <Rinternals.h>
#include <string.h> // for memset()
#include <stdint.h> // uint32_t, etc.

#define MAX_K 15

// Count all possible kmers for nucleic acid sequences
// Allow kmers up to 15 bases.
// Anything more and we might run into trouble
// since I will simply use an integer offset to represent kmers
//
// Although we can apparently return vectors longer than 2^31,
// if I want to deal with specific offsets in R I may run into trouble
// for example if I want to obtain the kmer associated with a given
// count.

// Note that bit operations are only defined on unsigned integers

// char -> offset is done by the following logic:
/* 
A  41  -> 0100 0001
C  43  -> 0100 0011
G  47  -> 0100 0111
T  54  -> 0101 0100

mask (0x06) everything except bits 2 and 3 and shift one to the right
Then simply shift the mask two steps to the left for each characther
and do a bitwise OR

Note that this does not check for Ns and other ambiguity symbols
We could include a check for this, but if we do not expect a very
large number of such letters then this may not matter. I should probably
make two functions. One that checks and one that doesn't bother.

N  4E  -> 0100 1110
this means that N will be counted as G (11).
*/

void increment_kmer_count(const char *seq, int *counts, uint8_t k){
  const char *s = seq;
  size_t r_mask = 0x06;  // right mask
  size_t l_mask = (1 << (k * 2)) - 1;  // left mask
  uint32_t offset = 0;
  uint8_t i = 0;
  for(i=0; i < k && *s; ++i){
      offset = (offset << 2) | ((*s & r_mask) >> 1);
      ++s;
  }
  if(i < k)
    return;
  ++counts[offset];
  // and then the remaining ones
  while(*s){
    offset = l_mask & ((offset << 2) | ((*s & r_mask) >> 1));
    ++counts[offset];
    ++s;
  }
}

// the memory for the kmer must be assigned before hand and must
// have length k+1, and ending in a 0.
// This assumes 
void int_to_kmer( uint32_t offset, uint8_t k, char *kmer ){
  char nuc[4] = {'A', 'C', 'T', 'G'};
  offset <<= ( 32 - (k * 2) );
  for(uint8_t i=0; i < k; ++i){
    kmer[i] = nuc[ offset >> 30 ];
    offset <<= 2;
  }
}

void print_bin(uint32_t value){
  uint32_t left_mask = 1 << 31;
  for(int i=0; i < 32; ++i){
    Rprintf("%d", (value & left_mask) >> 31);
    value <<= 1;
  }
  Rprintf("\n");
}

// offset XOR 0xAA 
// will complement. 
// Reverse bit pairs in k/2 operations..
// This seems to work, but could do with a cleanup and some better
// explanation of what the function does.
uint32_t rev_comp( uint32_t offset, uint8_t k ){
  uint32_t mask = (1 << k*2) - 1;
  offset = mask & (offset ^ 0xAAAAAAAA); // complements the sequence
  uint32_t nuc_mask = 3; // 0xb11  can be shifted appropriate number of times to mask unchanged positions.
  int8_t shift = k-1;
  while( shift >= k/2 ){
    uint32_t right = nuc_mask << (2 * (k - (shift+1)));
    uint32_t left = nuc_mask << (2 * shift);
    uint8_t dst = (2 * shift) - (2 * (k - (shift+1)));
    offset = (offset & ~(right | left)) | ((offset & right) << (dst)) | ((offset & left) >> (dst));
    shift -= 1;
  }
  offset = offset & mask;
  return(offset);
}

// A function to make the kmer counts from a character
// vector and a kmer size.
SEXP count_kmers( SEXP seqs_r, SEXP k_r ){
  if(TYPEOF(seqs_r) != STRSXP || length(seqs_r) < 1)
    error("first argument should be a non-null character vector");
  if(TYPEOF(k_r) != INTSXP || length(k_r) != 1)
    error("second argument (k) should be a single integer");
  int k_i = INTEGER(k_r)[0];
  if( k_i < 1 || k_i > MAX_K )
    error("invalid k. k should be 0 < k <= %d", MAX_K);
  uint8_t k = (uint8_t)k_i;
  // Assign a suitable vector for the return data:
  size_t k_count_size = (1 << (2 * k));
  SEXP kcounts_r = PROTECT( allocVector(INTSXP, k_count_size) );
  int *kcounts = INTEGER(kcounts_r);
  memset( kcounts, 0, sizeof(int) * k_count_size );
  // Go through each of the sequences and assign.
  int seq_n = length( seqs_r );
  for(int i=0; i < seq_n; ++i){
    increment_kmer_count( CHAR(STRING_ELT(seqs_r, i)), kcounts, k);
  }
  UNPROTECT(1);
  return( kcounts_r );
}

SEXP ints_to_kmers(SEXP offsets_r, SEXP k_r){
  if(TYPEOF(offsets_r) != INTSXP || length(offsets_r) < 1)
    error("offsets must be a vector of ints of positive length");
  if(TYPEOF(k_r) != INTSXP || length(k_r) != 1)
    error("second argument (k) should be a single integer");
  int n = length(offsets_r);
  int k = INTEGER(k_r)[0];
  if(k < 1 || k > MAX_K)
    error("Invalid k value");
  int *offsets = INTEGER(offsets_r);
  SEXP ret_value = PROTECT(allocVector( STRSXP, n ));
  char *kmer = malloc((k + 1) * sizeof(char));
  kmer[k] = 0;
  for(int i=0; i < n; ++i){
    int_to_kmer( (uint32_t)offsets[i], (uint8_t)k, kmer );
    SET_STRING_ELT(ret_value, i, mkChar( (const char*)kmer ));
  }
  free(kmer);
  UNPROTECT(1);
  return(ret_value);
}

SEXP rc_kmer_ints(SEXP offsets_r, SEXP k_r){
  if(TYPEOF(offsets_r) != INTSXP || length(offsets_r) < 1)
    error("offsets must be a vector of ints of positive length");
  if(TYPEOF(k_r) != INTSXP || length(k_r) != 1)
    error("second argument (k) should be a single integer");
  int n = length(offsets_r);
  int k = INTEGER(k_r)[0];
  if(k < 1 || k > MAX_K)
    error("Invalid k value");
  int *offsets = INTEGER(offsets_r);
  SEXP ret_value = PROTECT(allocVector( INTSXP, n ));
  int *rc = INTEGER(ret_value);
  memset(rc, 0, sizeof(int) * n);
  for(int i=0; i < n; ++i){
    rc[i] = rev_comp( (uint32_t)offsets[i], (uint8_t)k );
  }
  UNPROTECT(1);
  return(ret_value);
}
	  
  
  
