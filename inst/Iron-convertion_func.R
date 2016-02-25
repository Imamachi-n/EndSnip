## Alpine core function: core.R <https://github.com/mikelove/alpine>
mapTxToGenome <- function(exons){
    #Get strand information of each gene
    strand <- as.character(strand(exons)[1])
    #Check exon numbers
    stopifnot(all(exons$exon_rank == seq_along(exons)))
    #each base list
    bases <- if (strand == "+") {
        do.call(c, lapply(exons, function(exon) start(exon):end(exon)))
    } else if (strand == "-") {
        do.call(c, lapply(exons, function(exon) end(exon):start(exon)))
    }
    #ID - bases(genome position) - exon#
    data.frame(tx = seq_along(bases),
               genome = bases,
               exon_rank=rep(exons$exon_rank, width(exons)))
}

genomeToTx <- function(genome, map) map$tx[match(genome, map$genome)]

txToGenome <- function(tx, map) map$genome[match(tx, map$tx)]

txToExon <- function(tx, map) map$exon_rank[match(tx, map$tx)]

startLeft <- function(x) {
    first.plus <- as.logical(strand(first(x)) == "+")
    ifelse(first.plus, start(first(x)), start(last(x)))
}

endRight <- function(x) {
    first.plus <- as.logical(strand(first(x)) == "+")
    ifelse(first.plus, end(last(x)), end(first(x)))
}

#ga <- ga
#grl <- GRangesList(gene)
#fco <- fco
gaToReadsOnTx <- function(ga, grl, fco = NULL) {
    reads <- list()
    for (i in seq_along(grl)) {
        exons <- grl[[i]] #GRanges: exons
        strand <- as.character(strand(exons)[1])
        #Extract read ID (Compatible reads)
        read.idx <- if (is.null(fco)) {
            seq_along(ga)
        } else {
            queryHits(fco)[subjectHits(fco) == i]
        }
        #ID - bases(genome position) - exon#
        map <- mapTxToGenome(exons)
        
        #Define start/end site on transcript
        if (strand == "+") {
            start <- genomeToTx(startLeft(ga[read.idx]), map)
            end <- genomeToTx(endRight(ga[read.idx]), map)
        } else if (strand == "-") {
            start <- genomeToTx(endRight(ga[read.idx]), map)
            end <- genomeToTx(startLeft(ga[read.idx]), map)
        }
        #Validation(Remove wrong paired-reads)
        vaild <- start < end & !is.na(start) & !is.na(end)
        reads[[i]] <- IRanges(start[vaild], end[vaild])
    }
    
    names(reads) <- names(grl)
    return(reads)
}



