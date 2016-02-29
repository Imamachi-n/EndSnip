## Alpine function: vlmm.R <https://github.com/mikelove/alpine>

#Output k-mer nucleotide patterns
alphafun <- function(dna.letters, order) {
    if (order == 0) {
        return(dna.letters)
    } else {
        return(as.vector(t(outer(alphafun(dna.letters, order-1), x, paste0))))
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

fitVLMM <- function(seqs, gene.seqs) {
    
}


