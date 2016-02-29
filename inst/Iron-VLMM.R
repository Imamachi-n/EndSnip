## Alpine function: vlmm.R <https://github.com/mikelove/alpine>

#Output k-mer nucleotide patterns
alphafun <- function(dna.letters, order) {
    if (order == 0) {
        return(dna.letters)
    } else {
        return(as.vector(t(outer(alphafun(dna.letters, order-1), dna.letters, paste0))))
    }
}


#Observed nucleotide frequency
getPositionalKmerFregs <- function(seqs, dna.letters, order, pos, pc=1) {
    #Make k-mer nucleotide patterns
    alpha <- alphafun(dna.letters, order)
    
    #Count each nucleotide patterns
    out <- as.numeric(table(factor(substr(seqs, pos-order, pos), alpha)))
    names(out) <- alpha
    out <- out + pc #Add pseudo-count
    return(out / sum(out)) #Return the percentage of each nucleotide patterns
}


#Expected nucleotide frequency
getKmerFreqs <- function(seqs, dna.letters, order, pc=1) {
    #Make k-mer nucleotide patterns
    alpha <- alphafun(dna.letters, order)
    
    #The number of 5'/3'side 21nt sequences
    n <- sum(width(seqs)) - order*length(seqs)
    
    #A,T,G,C frequency
    if (order > 0) {
        out <- sapply(alpha, function(p) sum(vcountPattern(p, seqs)))
    } else {
        out <- colSums(letterFrequency(seqs, dna.letters))
    }
    
    #Count each nucleotide patterns
    stopifnot(sum(out) == n)
    out <- out + pc #Add pseudo-count
    return(out / sum(out)) #Return the percentage of each nucleotide patterns
}


#Observed / Expected nucleotide frequency
getPositionalObsOverExp <- function(seqs, gene.seqs, dna.letters, order, pos) {
    #Position
    npos <- length(pos)
    
    ## -- Observed nucleotide frequency(1-order/2-order) --
    res <- sapply(pos, function(i) getPositionalKmerFregs(seqs, dna.letters, order = order, pos = i))
    
    #'obs' is a dimensional array:
    # 1st dim: A,C,G,T the current base
    # 2nd dim: the 4(A,C,G,T)^order previous bases
    # 3rd dim: the position, which is a subset of the full VLMM order
    
    #Initialize 3-dimensional array
    obs <- array(0, dim = c(4, 4^order, npos), dimnames = list(dna.letters, alphafun(dna.letters, order - 1), seq_len(npos)))
    
    #Observed nucleotide frequency(1nt->2nt/2nt->3nt)
    for (i in 1:npos) {
        #Copy each di/tri-nucleotides frequency at each position
        obs[,,i] <- res[,i]
        
        #Each nucleotide frequency(e.g. AA/(AA+AC+AG+AT), AC/(AA+AC+AG+AT))
        obs[,,i] <- sweep(obs[,,i], 2, colSums(obs[,,i]), "/")
    }
    
    ## -- Expected nucleotide frequency(1-order/2-order) --
    res.gene <- getKmerFreqs(gene.seqs, dna.letters, order)
    alpha <- alphafun(dna.letters, order - 1)
    
    #'expect' has dims 1 and 2 from above
    expect <- array(0, dim = c(4, 4^order), dimnames = list(dna.letters, alphafun(dna.letters, order - 1)))
    
    for (p in alpha) {
        #Grep 1nt/2nt nucleotides patterns ([AA,CA,GA,TA],[AC,CC,GC,TC]...)
        prob <- res.gene[grep(paste0("^", p), names(res.gene))]
        
        #Calculate each nucleotide frequency
        expect[,p] <- prob / sum(prob)
    }
    
    #Result
    return(list(obs = obs, expect = expect))
}


#Calculate Observed/Expected nucleotide frequency (0,1,2-order)
#seqs <- fivep
#gene.seqs <- gene.seqs
fitVLMM <- function(seqs, gene.seqs) {
    #Fit a VLMM according to Roberts et al. (2011), doi:10.1186/gb-2011-12-3-r22
    dna.letters <- c("A", "C", "G", "T")
    
    #Parameter
    vlmm.order <- c(0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 0, 0)
    
    ## -- 0-order --
    #Initialize list
    order0 <- list()
    
    #Observed nucleotide frequency (0-order)
    order0$obs <- sapply(seq_along(vlmm.order), function(i) getPositionalKmerFregs(seqs, dna.letters, order = 0, pos = i))
    colnames(order0$obs) <- seq_along(vlmm.order)
    
    #Expected nucleotide frequency (0-order)
    order0$expect <- getKmerFreqs(gene.seqs, dna.letters, 0)
    
    ## -- 1-order --
    order <- 1
    pos1 <- which(vlmm.order >= order)
    order1 <- getPositionalObsOverExp(seqs, gene.seqs, dna.letters, order, pos1)
    
    ## -- 2-order --
    order <- 2
    pos2 <- which(vlmm.order >= order)
    order2 <- getPositionalObsOverExp(seqs, gene.seqs, dna.letters, order, pos2)
    
    #Result
    return(list(order0 = order0, order1 = order1, order2 = order2))
}

calcVLMMBias <- function(seqs, vlmm.model, short = FALSE, pseudocount = 1) {
    
}

addVLMMBias <- function(fragtypes, vlmm.fivep, vlmm.threep) {
    
}





