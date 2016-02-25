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

txToGenome <- function(tx, map) map$genome[match(tx, map$tx)]
