source('paf_functions.R')
dyn.load("parse_paf.so")

paf <- read.paf('test.paf')
cs <- parse.cs(paf)
cs.mut <- lapply( cs, paf.mut )
cs.mut.ch <- lapply( cs, paf.mut.ch )

## with the integral version we can do:
mut.count <- apply( do.call(rbind, cs.mut), 2, table )

intToUtf8( names(mut.count$r) ) ##
intToUtf8( names(mut.count$q) ) ## 
## that gives me lots of 0s; which we shouldn't have


head(i <- grep('\\+', paf$cs))

## example with +
paf$cs[i[1]] ## :397+t:140
paf$cs[i[3]] ## :84*ag:86+g:34*gt:28*ga:68*gt:4*ag:36+a:8*ga:5*ag:5*ga:78


tmp <- .Call("parse_diff_string", paf$cs)

intToUtf8( tmp[[1]][,1] )

