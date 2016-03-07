## Alpine function: fit_bias.R <https://github.com/mikelove/alpine>

#FPBP needed to downsample to a target fragment per kilobase
getFPBP <- function(genes, bamfile){
    gene.ranges <- unlist(range(genes)) #GRanges: st -> ed (including exon/intron)
    gene.lengths <- sum(width(genes))
    
    #Count mapped reads on each transcript
    res <- countBam(bamfile, param = ScanBamParam(which = gene.ranges))
    
    #Normalized counts
    out <- (res$records / 2) / gene.lengths
    names(out) <- names(genes)
    return(out)
}

#Count mapped reads on each transcript
alpineFlag <- function() scanBamFlag(isSecondaryAlignment = FALSE)
readGAlignAlpine <- function(bamfile, generange, readType){
    if (readType == "SE") {
        readGAlignments(bamfile, param = ScanBamParam(which = generange, flag = alpineFlag()))
    } else if (readType == "PE") {
        readGAlignmentPairs(bamfile, param = ScanBamParam(which = generange, flag = alpineFlag()))
    } else {
        stop("readType is wrong. Select SE or PE.")
    }
}

#
matchToDensity <- function(x, d) {
    #Devide each fragment based on the distribution of fragment length
    idx <- cut(x, c(-Inf, d$x, Inf))
    
    #Prepare dentity score
    pdf <- c(0, d$y)
    
    #
    pdf.x <- pdf[idx] + 1e-6
    
    stopifnot(all(pdf.x > 0))
    return(pdf.x)
}

#GLM fitting

#genes <- ebt.rep
#bamfile <- bamfile
#genome <- Hsapiens
#models <- models
#readType <- "SE"
#readlength <- 36
#zerotopos <- 2
#speedglm <- TRUE
#minsize <- 100
#maxsize <- 300

fitModelOverGenes <- function(genes, bamfile, genome,
                              models, readType, readlength, zerotopos = 2,
                              speedglm = TRUE, minsize = 100, maxsize = 300) {
    #Checking
    stopifnot(file.exists(bamfile))
    stopifnot(file.exists(paste0(as.character(bamfile),".bai")))
    stopifnot(is(genes, "GRangesList"))
    stopifnot(all(!is.na(sapply(models, function(x) x$formula))))
    stopifnot(is.numeric(readlength) & length(readlength) == 1)
    
    
    #Extract transcript sequence
    exon.dna <- getSeq(genome, genes) #DNAStringSetList: Exon level for each gene(transcript)
    gene.seqs <- as(lapply(exon.dna, unlist), "DNAStringSet") #DNAStringSet: Transcript level for each gene(transcript)
    
    #FPBP needed to downsample to a target fragment per kilobase
    fpbp <- getFPBP(genes, bamfile)
    
    target.fpbp <- 0.2 #FPBP Criteria
    
    #Fitting parameters
    fitpar.sub <- list()
    fitpar.sub[["coefs"]] <- list()
    
    fragtypes.sub.list <- list()
    #i <- 1 #TEST:
    for (i in seq_along(genes)) {
        gene.name <- names(genes)[i]
        gene <- genes[[gene.name]] #GRanges: gene infor
        l <- sum(width(gene)) #length
        
        #Add counts for sample and subset
        generange <- range(gene) #GRanges
        strand(generange) <- '*' #not necessary
        
        #Checking: chromosome number is the same as genome
        if (!as.character(seqnames(generange)) %in% seqlevels(BamFile(bamfile))) next
        
        #Mapped reads on each transcript (st -> ed (including exon/intron))
        suppressWarnings({
            ga <- readGAlignAlpine(bamfile, generange, readType) #GAlignments object
        })
        
        #Remove genes with Low coverage 
        if (length(ga) < 20) next
        
        ga <- keepSeqlevels(ga, as.character(seqnames(gene)[1]))
        
        #Calculate normalized read counts
        nfrags <- length(ga)
        this.fpbp <- nfrags / l
        
        #Downsampling if sequence reads is too large, ... # TODO: Need to change the timing of downsampling ?? (After findCompatibleOverlaps)
        if (this.fpbp > target.fpbp) {
            ga <- ga[sample(nfrags, round(nfrags * (target.fpbp / this.fpbp)), FALSE)]
        }
        
        #Compatible reads on each transcript (st -> ed (including only exon!!))
        fco <- findCompatibleOverlaps(ga, GRangesList(gene)) #Hits object
        
        #Transcript position of compatible reads
        reads <- gaToReadsOnTx(ga, GRangesList(gene), fco, readType, readlength)
        
        #Prepare dummy data(All possible fragment patterns)
        fragtypes <- list()
        if (readType == "SE") {
            fragtypes <- list(buildFragtypesFromExonsSE(genes[[gene.name]], genome = Hsapiens, readlength = readlength, 
                                                        npre = 8, npost = 12, gc = TRUE, gc.str = TRUE, vlmm = TRUE))
        } else if (readType == "PE") {
            fragtypes <- list(buildFragtypesFromExons(genes[[gene.name]], genome = Hsapiens,
                                                      readlength = readlength, minsize = 100, maxsize = 300))
        }

        #Checking
        #stopifnot(all(names(genes) %in% names(fragtypes)))
        #if (any(sapply(models, function(m) "vlmm" %in% m$offset))) {
        #    stopifnot("fivep" %in% colnames(fragtypes[[1]]))
        #}
        
        #Add count data to fragtypes(dummy data)
        fraglist.temp <- matchReadsToFraglist(reads, fragtypes)
        
        #Remove the following reads: Start/End site == 1
        not.first.or.last.bp <- !(fraglist.temp[[1]]$start == 1 | fraglist.temp[[1]]$end == 1)
        fraglist.temp[[1]] <- fraglist.temp[[1]][not.first.or.last.bp,]
        
        #Check read depth
        if (sum(fraglist.temp[[1]]$count) < 20) next
        
        #Subset to include all (N) fragment locations and (N * zerotopos) zero locations
        #Now we build a list of subsetted fragtypes
        fragtypes.sub.list[[gene.name]] <- subsetAndWeightFraglist(fraglist.temp, zerotopos)
        
    }
    
    #Checking
    if (length(fragtypes.sub.list) == 0) stop("not enough reads to model: ", bamfile)
    if (length(fragtypes.sub.list) < 2) stop("requires at least two genes to fit model")
    
    #Check the number of fragments(mapped reads + dummy reads)
    gene.nrows <- sapply(fragtypes.sub.list, nrow)
    
    #Merge each gene data
    fragtypes.sub <- do.call(rbind, fragtypes.sub.list)
    
    #Prepare fitting parameters
    fitpar.sub[["models"]] <- models
    
    ## -- Fragment bias --
    if (any(sapply(models, function(m) "fraglen" %in% m$offset))) {
        #Select mapped reads from total reads(mapped reads/dummy reads)
        pos.count <- fragtypes.sub$count > 0
        
        #Selected fragment length from mapped reads    
        fraglens <- rep(fragtypes.sub$fraglen[pos.count], fragtypes.sub$count[pos.count])
        
        #Calculate the density of fragment length from mapped reads
        fraglen.density <- density(fraglens)
        
        fragtypes.sub$logdfraglen <- log(matchToDensity(fragtypes.sub$fraglen, fraglen.density))
        
        fitpar.sub[["fraglen.density"]] <- fraglen.density
    }
    
    ## -- Random hexamer priming bias with VLMM --
    if (sapply(models, function(m) "vlmm" %in% m$offset)) {
        if (readType == "SE") {
            #5'side sequence
            fivep <- fragtypes.sub$fivep[fragtypes.sub$fivep.test]
            
            #Observed/Expected nucleotide frequency (5'/3'side)
            vlmm.fivep <- fitVLMM(fivep, gene.seqs)
        } else if (readType == "PE") {
            #5'side sequence
            fivep <- fragtypes.sub$fivep[fragtypes.sub$fivep.test]
            #3'side sequence
            threep <- fragtypes.sub$threep[fragtypes.sub$threep.test]
            
            #Observed/Expected nucleotide frequency (5'/3'side)
            vlmm.fivep <- fitVLMM(fivep, gene.seqs)
            vlmm.threep <- fitVLMM(threep, gene.seqs)
        }

        #Now calculate log(bias) for each fragment based on the VLMM
        if (readType == "SE") {
            fragtypes.sub <- addVLMMBiasSE(fragtypes.sub, vlmm.fivep)
        } else if (readType == "PE") {
            fragtypes.sub <- addVLMMBias(fragtypes.sub, vlmm.fivep, vlmm.threep)
        }

        #Add VLMM parameters into 'fitpar.sub'
        fitpar.sub[["vlmm.fivep"]] <- vlmm.fivep
    }
    
    #Allow a gene-specific intercept (although mostly handled already with downsampling)
    #Add geneID into 'fragtypes.sub'
    fragtypes.sub$gene <- factor(rep(seq_along(gene.nrows), gene.nrows))
    
    #modeltype <- "GC"
    for (modeltype in names(models)) {
        #Parameters
        #GC-contents
        gc.knots <- seq(from = .4, to = .6, length = 3) #[0.4, 0.5, 0.6]
        gc.bk <- c(0, 1) #[0, 1]
        #The Position of fragments
        relpos.knots <- seq(from = .25, to = .75, length = 3) #[0.25, 0.50, 0.75]
        relpos.bk <- c(0, 1) #[0, 1]
        
        #Prepare formula
        f <- models[[modeltype]]$formula
        
        #Initialize offset
        offset <- numeric(nrow(fragtypes.sub))
        
        ## -- Fragment bias --
        if ("fraglen" %in% models[[modeltype]]$offset) {
            offset <- offset + fragtypes.sub$logdfraglen
        }
        
        ## -- Random hexamer priming bias with VLMM --
        if ("vlmm" %in% models[[modeltype]]$offset) {
            if (readType == "SE") {
                offset <- offset + fragtypes.sub$fivep.bias
            } else if (readType == "PE") {
                offset <- offset + fragtypes.sub$fivep.bias + fragtypes.sub$threep.bias
            }
        }
        
        #Checking
        if (!all(is.finite(offset))) stop("offset needs to be finite")
        
        #Add offset into 'fragtypes.sub'
        fragtypes.sub$offset <- offset
        
        if (speedglm) {
            #Design matrix
            mm.small <- model.matrix(formula(f), data = fragtypes.sub)
            
            #Checking
            stopifnot(all(colSums(abs(mm.small) > 0)))
            
            #GLM fitting(speedglm)
            fit <- speedglm.wfit(fragtypes.sub$count, mm.small,
                                 family = poisson(),
                                 weights = fragtypes.sub$wts,
                                 offset = fragtypes.sub$offset)
        } else {
            fit <- glm(formula(f), family = poisson, data = fragtypes.sub, weights = wts, offset = offset)
        }
        fitpar.sub[["coefs"]][[modeltype]] <- fit$coefficients
        fitpar.sub[["summary"]][[modeltype]] <- summary(fit)$coefficients
    }
    return(fitpar.sub)
}