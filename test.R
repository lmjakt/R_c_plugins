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


### try the watershed segmentation
dyn.load("water_shed_segment.so")

mat <- readRDS("blurred_matrix.rds")
cols <- hcl.colors(256, "YlOrRd", rev=TRUE)
mat.t <- mat
mat.t <- mat.t - min(mat.t)
x <- 1:nrow(mat) - 1
y <- 1:ncol(mat) - 1
par(mfrow=c(1,3))
image(x, y, mat, col=cols)
##image(mat.t, col=cols)
mat.seg <- .Call("ws_segment", mat, quantile(mat)[4])
points(mat.seg[[2]][,1], mat.seg[[2]][,2], col='light blue')
image(x, y, mat.seg[[1]])
points(mat.seg[[2]][,1], mat.seg[[2]][,2], col='blue')
## first column is the row, but this gets turned to
## be the columns..

ridges <- .Call("ws_trace_ridges", mat, mat.seg[[2]], 0.0)
image(x, y, ridges[[1]])
points(mat.seg[[2]][,1], mat.seg[[2]][,2], col='blue', pch=19)

hist2d <- readRDS("hist2d.rds")
with(hist2d, image(x, y, h), col=cols)


## this meas that we can classify all the points by the mat.seg..
pts.id <- with( hist2d, sapply(1:length(xv), function(i){
    mat.seg[[1]][ 1 + xv[i] - x[1], 1 + yv[i] - y[1] ]
}))

par(mfrow=c(1,3))
with(hist2d, plot(xv, yv, col=rgb(0,0,0,0.2)))
with(hist2d, image(x, y, mat.seg[[1]]))
with(hist2d, {b <- pts.id == 1; points( xv[b], yv[b], cex=0.4 ) })
with(hist2d, {b <- pts.id == 2; points( xv[b], yv[b], cex=0.4, col='blue' ) })
with(hist2d, {b <- pts.id == 3; points( xv[b], yv[b], cex=0.4, col='grey' ) })
with(hist2d, {b <- pts.id == 4; points( xv[b], yv[b], cex=0.4, col='green' ) })
with(hist2d, {b <- pts.id == 5; points( xv[b], yv[b], cex=0.4, col='yellow' ) })

with(hist2d, {b <- pts.id == 1; lines(lowess( xv[b], yv[b]), col='green', lwd=3)})
with(hist2d, {b <- pts.id > 1; lines(lowess( xv[b], yv[b]), col='green', lwd=3)})
## lv <- with(hist2d, {b <- pts.id == 3; lines(lowess( xv[b], yv[b]), col='green', lwd=3)})
## lv <- with(hist2d, {b <- pts.id == 4; lines(lowess( xv[b], yv[b]), col='green', lwd=3)})
## lv <- with(hist2d, {b <- pts.id == 5; lines(lowess( xv[b], yv[b]), col='green', lwd=3)})

with(hist2d, image(x, y, mat, col=cols))

require('entropy')
## just to see how discretize2d works (uses cut; a very useful function)

### deivide the range of x into intervals..
x.r <- range(hist2d$x)
x.b <- seq(from=x.r[1], to=x.r[2], length.out=11)
x.b2 <- x.b + diff(x.b)[1] / 2
x.l <- cut(hist2d$xv, x.b, include.lowest=TRUE)
x.l2 <- cut(hist2d$xv, x.b2, include.lowest=TRUE)
y.d <- tapply( hist2d$yv, x.l, density )
y.d2 <- tapply( hist2d$yv, x.l2, density )

i.col <- hcl.colors(length(y.d), "YlOrRd")
plot(y.d[[1]], col=i.col[[10]])
with(par(), rect(usr[1], usr[3], usr[2], usr[4], col='grey'))
for(i in 1:length(y.d)){
    lines(y.d[[i]], col=i.col[i])
    lines(y.d2[[i]], col=i.col[i])
}

with(y.d[[1]], plot(x[-1], diff(y)))
abline(h=0, lty=2)
plot(y.d[[1]])

for(k in 1:length(y.d)){
    ## this is ugly, but..
    i <- with(y.d[[k]], {
        d <- diff(y)
        which(sapply(4:(length(d)-3), function(i){
            all(c( d[(i-3):(i-1)] > 0, d[(i+1):(i+3)] < 0))
        }))
    })
    ##
    par(mfrow=c(1,1))
    plot(y.d[[k]], type='b', cex=0.5) ## , xlim=c(50, 70))
    abline(v=y.d[[k]]$x[ i + 4 ], lty=2)
    inpt <- readline("next: ")
    i <- with(y.d2[[k]], {
        d <- diff(y)
        which(sapply(4:(length(d)-3), function(i){
            all(c( d[(i-3):(i-1)] > 0, d[(i+1):(i+3)] < 0))
        }))
    })
    ##
    par(mfrow=c(1,1))
    plot(y.d2[[k]], type='b', cex=0.5) ## , xlim=c(50, 70))
    abline(v=y.d2[[k]]$x[ i + 4 ], lty=2)
    inpt <- readline("next: ")
}


with(y.d[[1]], plot(x[-1], diff(y), cex=0.5))
abline(h=0, lty=2)
abline(v=y.d[[1]]$x[ i + 4 ], lty=2)

get.peaks <- function(x, y, w){
    d <- diff(y)
    p.i <- which(sapply((w+1):(length(d)-w), function(i){
            all(c( d[(i-w):(i-1)] > 0, d[(i+1):(i+w)] < 0))
    }))
    if(length(p.i) < 2)
        return(p.i)
    ## we expect each peak to be represented by a doublet.
    p.d <- diff(p.i)
    p.db <- p.d == 1
    peak.pos <- c()
    i <- 1
    while(i <= length(p.i)){
        if(i < length(p.i) && p.i[i+1] - p.i[i] == 1){
            peak.pos <- c(peak.pos, (x[p.i[i]+w] + x[p.i[i+1]+w]) / 2)
            i <- i + 2
            next
        }
        peak.pos <- c(peak.pos, x[p.i[i]+w])
        i <- i + 1
    }
    peak.pos
}
        
    
    

## width and step are in terms of number of measures
## which are ordered by x
windowed.density <- function(x, y, width, step, max.d=0.02 * diff(range(y))){
    o <- order(x);
    beg <- seq(1, length(x)-width, step)
    end <- beg + width
    wd <- lapply( 1:length(beg), function(i){
        j <- o[beg[i]:end[i]]
        x.q <- quantile(x[j])
        x.m <- mean(x[j])
        d <- density( y[j] )
        p <- get.peaks( d$x, d$y, 3 )
        list(i=j, x.q=x.q, x.m=x.m, d=d, p=p)
    })
    ## merge peaks into lines..
    peaks <- lapply( wd, function(x){ cbind(x=x$x.m, p=x$p, assigned=FALSE) })
    ## a little bit of recursion might work here..
    trace.ridge <- function(i, j, ridge, max.d){
        if(i > length(peaks) || peaks[[i]][j,'assigned'])
            return(ridge)
        ridge <- rbind(ridge, peaks[[i]][j,1:2])
        peaks[[i]][j,'assigned'] <<- TRUE
        if(i == length(peaks))
            return(ridge)
        d <- abs( peaks[[i]][j,'p'] - peaks[[i+1]][,'p'] )
        d.i <- which.min(d)
        if(d[d.i] < max.d && !peaks[[i+1]][d.i,'assigned'] )
            ridge <- trace.ridge(i+1, d.i, ridge, max.d)
        ridge
    }
    ridges <- list()
    for(i in 1:length(peaks)){
        for(j in 1:nrow(peaks[[i]])){
            if(!peaks[[i]][j,'assigned']){
                ridge <- matrix(nrow=0, ncol=2)
                ridge <- trace.ridge(i, j, ridge, max.d)
                ridges <- c(ridges, list(ridge))
            }
        }
    }
    list(wd=wd, ridges=ridges)
}

dyn.load("gs_smooth.so")
    
tmp <- windowed.density( hist2d$xv, hist2d$yv, 250, 50, max.d=0.1*diff(range(hist2d$y)) )

min.f <- 0.01
sd <- 3.0
ridge.sm <- lapply( tmp$ridges, function(r){
    .Call("gs_smooth", r, sd, min.f)
})

with(hist2d, image(x, y, mat, col=cols))
lines( tmp$ridges[[1]][,1], ridge.sm[[1]], lwd=3, col='black' )
lines( tmp$ridges[[2]][,1], ridge.sm[[2]], lwd=2, col='black' )

lines( tmp$ridges[[2]][,1], ridge.sm[[2]], lwd=2, col='black' )


with(hist2d, image(x, y, h))
lines( tmp$ridges[[1]][,1], ridge.sm[[1]], lwd=3, col='green' )
ines( tmp$ridges[[2]][,1], ridge.sm[[2]], lwd=3, col='blue' )

points( tmp$ridges[[32]][,1], ridge.sm[[32]], lwd=3, col='blue' )
points( tmp$ridges[[32]][,1], ridge.sm[[32]], lwd=3, col='blue' )


with(hist2d, image(x, y, h))
lines(tmp$ridges[[1]], type='b')
lines(tmp$ridges[[2]], type='b')
