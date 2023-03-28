
## this performs some magic that should not be necessary if created
## as a package
dyn.load( paste(dirname(sys.frame(1)$ofile), "parse_paf.so", sep="/") )

require('parallel')

## this can be done much more elegantly, but ..
## this will complement
comp.table <- 1:128
comp.table[ utf8ToInt('A') ] <- utf8ToInt('T')
comp.table[ utf8ToInt('C') ] <- utf8ToInt('G')
comp.table[ utf8ToInt('G') ] <- utf8ToInt('C')
comp.table[ utf8ToInt('T') ] <- utf8ToInt('A')
##
comp.table[ utf8ToInt('a') ] <- utf8ToInt('t')
comp.table[ utf8ToInt('c') ] <- utf8ToInt('g')
comp.table[ utf8ToInt('g') ] <- utf8ToInt('c')
comp.table[ utf8ToInt('t') ] <- utf8ToInt('a')


read.paf <- function(fn){
    paf.lines <- strsplit( readLines(fn), "\t", fixed=TRUE)
    paf <- as.data.frame( t(sapply( paf.lines, function(x){ x[1:12] })), stringsAsFactors=FALSE )
    colnames(paf) <- c('query', 'qlen', 'qstart', 'qend', 'strand', 'ref', 'rlen', 'rstart', 'rend', 'match', 'al.l', 'mqual')
    col.num <- c(2:4, 7:12)
    for(i in col.num)
        paf[,i] <- as.numeric(paf[,i])
    ## obtain the values of de, cs, AS, ms, NM, nn if present
    tag.v <- sapply( paf.lines, function(x){
        cols <- sapply( strsplit( x[13:length(x)], "[AifZ]:", fixed=FALSE ), eval)
        if(is.null(dim(cols))){
            print(cols)
            stop("Dim of cols is wrong")
        }
        AS <- (cols[ 2, cols[1,] == 'AS:' ])[1]
        ms <-  (cols[ 2, cols[1,] == 'ms:' ])[1]
        de <-  (cols[ 2, cols[1,] == 'de:' ])[1]
        NM <-  (cols[ 2, cols[1,] == 'NM:' ])[1]
        nn <-  (cols[ 2, cols[1,] == 'nn:' ])[1]
        tp <- (cols[ 2, cols[1,] == 'tp:' ])[1]
        cs <- (cols[ 2, cols[1,] == 'cs:' ])[1]
        c(AS=AS, ms=ms, de=de, NM=NM, nn=nn, tp=tp, cs=cs)
    })
    tag.df <- as.data.frame( t(tag.v ), stringsAsFactors=FALSE )
    for(i in 1:5)
        tag.df[,i] <- as.numeric(tag.df[,i])
    cbind(paf, tag.df, stringsAsFactors=FALSE)
}

## takes an object returned by the read.paf function
## requires that dyn.load() has been called on parse_paf.so
## I should really do a R_registeroutines call in the data
## set.
parse.cs <- function(paf, include.strings=FALSE){
    tmp <- .Call("parse_diff_string", paf$cs, as.integer(paf$rstart),
                 as.integer(paf$qstart), as.integer(paf$qend),
                 paf$strand == '+', include.strings);
    names(tmp) <- c("al", "seq")
    if(!include.strings)
        return(tmp[[1]])
    tmp
}

## makes a data.frame with the operations and mutations shown..
cs.df <- function(cs, n=nrow(cs)){
    mut <- paf.mut(cs)
    mut.c <- apply(mut, 2, function(x){ strsplit(intToUtf8(x), "")[[1]] })
    ops <- strsplit( intToUtf8(cs[,'op']), "")[[1]]
    data.frame( opc=ops, mutc.r=mut.c[,'r'], mutc.q=mut.c[,'q'], cs, stringsAsFactors=FALSE)
}

## extract the mutations from the mut column
## takes a single entry from the list returned by
## parse.cs
## to integer values (more efficient)
paf.mut <- function(tbl){
    cbind('r'=bitwShiftR(tbl[,'mut'], 8),
          'q'= bitwAnd(tbl[,'mut'], 0xFF) )
}

## to character vectors; less efficient..
paf.mut.ch <- function(tbl){
    cbind('r'=intToUtf8(bitwShiftR(tbl[,'mut'], 8)),
          'q'=intToUtf8(bitwAnd(tbl[,'mut'], 0xFF)) )
}

paf.op <- function(tbl){
    intToUtf8(tbl[,'op'])
}


## takes a list of table and count the total
## set of mutations listed
paf.mut.count <- function(cs, strand, include.dot=FALSE, rev.comp=TRUE){
    tbl <- do.call(rbind, lapply(cs, paf.mut))
    if(rev.comp){
        rc <- ifelse(strand == '-', TRUE, FALSE)
        tbl[rc,1] <- comp.table[ tbl[rc,1] ]
        tbl[rc,2] <- comp.table[ tbl[rc,2] ]
    }
    counts <- table(tbl[,1], tbl[,2])
    colnames(counts) <- strsplit( intToUtf8(colnames(counts)), "")[[1]]
    rownames(counts) <- strsplit( intToUtf8(rownames(counts)), "")[[1]]
    if(include.dot)
        return(counts)
    counts[-1,-1]
}

make.paf.index <- function(paf, w.size=1000, primary.only=TRUE, min.ref=10000, cores=10){
    if(primary.only)
        i <- which( paf$tp == 'P' & paf$rlen >= min.ref)
    else
        i <- which( paf$rlen >= min.ref )
    ind <- tapply(i, paf$ref[i], eval)
    index <- mclapply( ind, function(j){
        chr.l <- paf$rlen[j[1]]
        start <- paf$rstart[j]
        end <- paf$rend[j]
        w.beg <- seq(0, chr.l, w.size)
        al.wi <- do.call(rbind, lapply(j, function(k){
            w.i <- as.integer(paf$rstart[k]/w.size):as.integer(paf$rend[k]/w.size)
            cbind(k, 1 + w.i)
        }))
        al.wi <- rbind( cbind(0, 1:length(w.beg)), al.wi )
        ll <- list('ind'=tapply(al.wi[,1], al.wi[,2], eval),
             'beg'=w.beg, 'end'=w.beg+w.size)
        ll[ order(as.numeric(names(ll))) ]
    }, mc.cores=cores)
    list(ind=index, w.size=w.size)
}

## takes an index, a reference and a position
## gives the index of all alignments that either include
## the relevant window, or have a start or end within in
## this means that some of the alingments may not include
## the position given... 
aligned.at <- function(ind, chr, pos){
    rows <- unique( ind$ind[[chr]]$ind[[ 1 + pos / ind$w.size ]] )
    rows[ rows != 0 ]
}

