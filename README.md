# A random selection of R plugins

With very little documentation.

## convolve_ks.c

Provides a naive symmetric kernel convolution function, `convolve_ks` that can be used
for various types of blurring and smoothing. Although implemented naively this function
appears to be much faster than the built-in `convolve` function. I have implemented
too many versions of this function and then forgotten where I put them; hope this is
the last one.

## nuc_kmer_count.c

Count nucleotide kmers in DNA or RNA strings. Fast function that converts k-mers to
unsigned integers; supports counting kmers up to 16 bases. Performs no error-checking and
returns a vector of counts for all possible kmers.

## spread_points.c

A function to spread points in an orthogonal dimension to decrease crowding of points.
Useful for 1 dimensional scatter plots.

## Usage

For each source file <source.c>:

``` sh
R CMD SHLIB <source.c>
```

And then `dyn.load` the resulting `source.so`. Read the `test.R` and the `c` sources
to work out how to use the functions.

