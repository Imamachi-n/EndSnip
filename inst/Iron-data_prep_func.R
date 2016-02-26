#exons <- ebt[[genenames[1]]]
#genome <- Hsapiens
#readlength <- 48
#minsize <- 100
#maxsize <- 300
#npre <- 8
#npost <- 12
#gc <- TRUE
#gc.str <- TRUE
#vlmm <- TRUE

buildFragtypesFromExons <- function(exons, genome, readlength,
                                    minsize, maxsize, npre = 8, npost = 12,
                                    gc = TRUE, gc.str = TRUE, vlmm = TRUE){
    #Checking
    stopifnot(is(exons, "GRanges"))
    stopifnot(is(genome, "BSgenome"))
    stopifnot(is.numeric(minsize) & is.numeric(maxsize) & is.numeric(readlength))
    stopifnot(sum(width(exons)) > maxsize)
    stopifnot(all(c("exon_rank", "exon_id") %in% names(mcols(exons))))
    
    #ID - bases(genome position) - exon#
    map <- mapTxToGenome(exons)
    
    #length of transcript
    l <- nrow(map)
    strand <- as.character(strand(exons)[1])
    
    #Fragment sites on transcripts
    start <- rep(seq_len(l - minsize + 1), each = maxsize - minsize + 1)
    end <- as.integer(start + minsize:maxsize - 1)
    
    #Relative positions on transcripts
    mid <- as.integer(0.5 * (start + end))
    rel_pos <- mid / l
    
    #Fragment length
    fraglen <- as.integer(end - start + 1)
    
    #IRanges: st - ed - width
    id <- IRanges(start, end)
    
    fragtypes <- DataFrame(start = start, end = end, rel_pos = rel_pos, fraglen = fraglen, id = id)
    #Filtering: unreliable fragments
    fragtypes <- fragtypes[fragtypes$end <= l,, drop=FALSE]
    
    #Get DNA sequences from genome
    exon.dna <- getSeq(genome, exons) #Exon level
    tx.dna <- unlist(exon.dna) #Transcript level
    
    ###Variable length Markov model (VLMM) for the random hexamer priming bias
    if (vlmm) {
        #5'side - sequenced read
        fragtypes$fivep.test <- fragtypes$start - npre >= 1 #Logical: From 9nt => extention
        fragtypes$fivep <- as(Views(tx.dna,
                                    fragtypes$start - ifelse(fragtypes$fivep.test, npre, 0),
                                    fragtypes$start + npost),
                              "DNAStringSet"
        ) #13nt or 21nt
        
        #3'side - sequenced read
        fragtypes$threep.test <- fragtypes$end + npost <= length(tx.dna) #'l' is also OK.
        fragtypes$threep <- as(Views(tx.dna,
                                     fragtypes$end - npost,
                                     fragtypes$end + ifelse(fragtypes$threep.test, npre, 0)),
                               "DNAStringSet"
        )
        
        #Reverse complement the three prime sequence
        fragtypes$threep <- reverseComplement(fragtypes$threep)
    }
    
    ###PCR amplification bias (Related with GC-contens)
    if (gc) {
        fragrange <- minsize:maxsize
        #Calcuate GC-contens for each fragment
        gc.vecs <- lapply(fragrange, function(i) {
            letterFrequencyInSlidingView(tx.dna, view.width = i, letters = "CG", as.prob = TRUE)
        })
        #Sort - fragment length to merge gc.vecs
        fragtypes <- fragtypes[order(fragtypes$fraglen),, drop=FALSE]
        fragtypes$gc <- do.call(c, gc.vecs)
        
        #Sort - fragment start site
        fragtypes <- fragtypes[order(fragtypes$start),, drop=FALSE]
    }
    
    ###Additional features: GC-contents in smaller sections
    #Window size: 40nt
    gc.40 <- as.numeric(letterFrequencyInSlidingView(tx.dna, 40, letters = "CG", as.prob = TRUE))
    #Select 5'/3' side sequence read with larger GC-contents
    max.gc.40 <- max(Views(gc.40,
                           fragtypes$start,
                           fragtypes$end - 40 + 1))
    
    #Window size: 20nt
    gc.20 <- as.numeric(letterFrequencyInSlidingView(tx.dna, 20, letters = "CG", as.prob = TRUE))
    #Select 5'/3' side sequence read with larger GC-contents
    max.gc.20 <- max(Views(gc.20,
                           fragtypes$start,
                           fragtypes$end - 20 + 1))
    
    #Result
    fragtypes$GC40.90 <- as.numeric(max.gc.40 >= 36/40) # 40nt - 90%
    fragtypes$GC40.80 <- as.numeric(max.gc.40 >= 32/40) # 40nt - 80%
    fragtypes$GC20.90 <- as.numeric(max.gc.20 >= 18/20) # 20nt - 90%
    fragtypes$GC20.80 <- as.numeric(max.gc.20 >= 16/20) # 20nt - 80%
    
    #Covert transcript to genome position
    fragtypes$gstart <- txToGenome(fragtypes$start, map)
    fragtypes$gend <- txToGenome(fragtypes$end, map)
    fragtypes$gread1end <- txToGenome(fragtypes$start + readlength - 1, map)
    fragtypes$gread2start <- txToGenome(fragtypes$end - readlength + 1, map)
    
    return(fragtypes)
}


#TEST: matchReadsToFraglist
#reads <- reads
#fraglist <- fragtypes[gene.name]
#tx.idx <- 1 #TEST:

matchReadsToFraglist <- function(reads, fraglist) {
    for (tx.idx in seq_along(fraglist)) {
        #Extract unique reads
        unique.reads <- unique(reads[[tx.idx]])
        
        #Count the number of reads at each position
        readtab <- table(match(reads[[tx.idx]], unique.reads))
        
        #Add 'count' column in "fraglist" DataFrame
        fraglist[[tx.idx]]$count <- 0
        
        #Search the same reads from dummy data(fragment range: 100-300nt)
        reads.in.fraglist <- unique.reads %in% fraglist[[tx.idx]]$id
        
        #Select unique reads matched with criteria(fragment range: 100-300nt)
        unique.reads <- unique.reads[reads.in.fraglist]
        readtab <- readtab[reads.in.fraglist]
        
        fraglist[[tx.idx]][match(unique.reads, fraglist[[tx.idx]]$id), "count"] <- as.numeric(readtab)
    }
    return(fraglist)
}


#TEST: subsetAndWeightFraglist
#fraglist <- fraglist.temp
#zerotopos <- 2
#minzero <- 2000
#maxmult <- 20000
#tx.idx <- 1 #TEST:

subsetAndWeightFraglist <- function(fraglist, zerotopos, minzero = 2000, maxmult = 20000) {
    ntx <- length(fraglist)
    fraglist.sub <- list()
    for (tx.idx in seq_len(ntx)) {
        #Total read counts
        count <- fraglist[[tx.idx]]$count
        
        sumpos <- sum(count > 0) #Mapped reads
        sumzero <- sum(count == 0) #Un-mapped reads
        
        #How many zeros of subset?
        #Some multiple of the #pos, or a preset minimum value
        multzero <- max(round(zerotopos * sumpos), minzero) #>=2000
        
        #Max out at a certain amount, and then sample equal to #pos
        multzero <- if (multzero > maxmult) max(maxmult, sumpos) else multzero
        
        #and not more than the #zero
        numzero <- min(sumzero, multzero)
        
        #Get selected ID after downsampling
        idx <- c(which(count > 0), 
                 sample(which(count == 0), numzero, replace = FALSE))
        
        #Subset the rows
        fraglist.sub[[tx.idx]] <- fraglist[[tx.idx]][idx,]
        
        #Add weights(Mapped reads->1 / Unmapped reads->SumZero/NumZero(normalized))
        fraglist.sub[[tx.idx]]$wts <- c(rep(1, sumpos), rep(sumzero/numzero, numzero))
        
        fragtypes <- do.call(rbind, fraglist.sub)
        return(fragtypes)
    }
}


