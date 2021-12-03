## kmer counting
seqs <- c("ACATACATACATACATACATGGGG",
          "ATAGACATAGAAGATAACATAGACGAGAACGAGA",
          "GGCGGCGGCGGCGGCGGCGGCGGCGGCGGC",
          "ccccccccccccccccccccccccccccccccccccccccc")
k <- 5L
dyn.load("nuc_kmer_count.so")
counts <- .Call("count_kmers", seqs, k)
kmers <- .Call("ints_to_kmers", as.integer( (1:length(counts))-1 ), k)
kmer.df <- data.frame(kmers, counts)
kmer.df <- kmer.df[ order(kmer.df$counts, decreasing=TRUE), ]

rev.comp <- .Call("rc_kmer_ints", as.integer( (1:length(counts))-1 ), k)

tmp <- .Call("rc_kmer_ints", as.integer( 21845 ), k); ## CCCCCCCC
tmp <- .Call("rc_kmer_ints", as.integer( 21846 ), k); ## CCCCCCCT
.Call("ints_to_kmers", as.integer(tmp), k )

## kernel convolution
dyn.load("convolve_ks.so")

values = as.double(sample(c(0,100), size=50, prob=c(0.8, 0.2), replace=TRUE))
## a square kernel

ks <- as.double(rep(1, 5))
ks.c <- .Call( "convolve_ks", values, ks )
plot(values)
lines(1:length(ks.c), ks.c)


ks <- dnorm(-2:2)
ks.c <- .Call( "convolve_ks", values, ks )
plot(values)
lines(1:length(ks.c), ks.c)


## a bigger data set:

values = as.double(sample(c(0,100), size=10000, prob=c(0.98, 0.02), replace=TRUE))
ks <- dnorm(0:40, sd=8)
ks.c <- .Call( "convolve_ks", values, ks )
plot(values)
lines(1:length(ks.c), ks.c * 10)

## and with a bigger terminal:
values = as.double(sample(c(0,100), size=50000, prob=c(0.995, 0.005), replace=TRUE))
ks <- dnorm(0:60, sd=20)
ks.c <- .Call( "convolve_ks", values, ks )
plot(values, xlim=c(1,1000))
lines(1:length(ks.c), ks.c * 30)


