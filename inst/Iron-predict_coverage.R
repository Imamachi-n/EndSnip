## Alpine function: predict.R <https://github.com/mikelove/alpine>

#gene <- ebt2[[test_ELAVL1_RXID]]
#genome <- Hsapiens
#readlength <- 48
#minsize <- 100
#maxsize <- 300

predictOneGene <- function(gene, bamfile, fitpar, genome=Hsapiens,
                           models, readlength, minsize, maxsize) {
    #Checking
    stopifnot(is(gene, "GRanges"))
    stopifnot(all(sapply(models, function(x) names(x) %in% c("formula","offset"))))
    stopifnot(!is.null(fitpar))
    
    #bamfile <- list.files(bamfile)
    
    #Prepare dummy reads(All possible fragment patterns)
    fragtypes <- buildFragtypesFromExons(gene, genome, readlength=readlength,
                                         minsize=minsize, maxsize=maxsize)
    
    #Initialize result list
    res <- list()
    
    #The range of test_gene (GRanges object)
    generange <- range(gene)
    strand(generange) <- "*" # not necessary
    
    #Mapped reads on each transcript (st -> ed (including exon/intron))
    suppressWarnings({
        ga <- readGAlignAlpine(bamfile, generange)
    })
    
    #IF: mapped reads -> 0
    if (length(ga) == 0) {
        res[["test"]] <- as.list(rep(NA,length(models)))
        names(res[["test"]]) <- names(models)
        next
    }
    
    ga <- keepSeqlevels(ga, as.character(seqnames(gene)[1]))
    
    #Compatible reads on each transcript (st -> ed (including only exon!!))
    fco <- findCompatibleOverlaps(ga, GRangesList(gene))
    
    #Transcript position of compatible reads
    reads <- gaToReadsOnTx(ga, GRangesList(gene), fco)
    
    # save fragment coverage for later
    l <- sum(width(gene)) #Transcript length
    
    #Remove the following reads: Start/End site == 1
    frag.cov <- coverage(reads[[1]][start(reads[[1]]) != 1 & end(reads[[1]]) != l])
    
    #Add count data to fragtypes(dummy data)
    fragtypes.temp <- matchReadsToFraglist(reads, list(fragtypes))[[1]]
    
    modeltype <- "GC"
    
    ## -- fragment bias --
    if ("fraglen" %in% models[[modeltype]]$offset) {
        fraglen.density <- fitpar[[1]][["fraglen.density"]]
        stopifnot(!is.null(fraglen.density)) #Checking
        
        #Calculate the density of fragment length from mapped reads
        fragtypes.temp$logdfraglen <- log(matchToDensity(fragtypes.temp$fraglen, fraglen.density))
    }

    ## -- random hexamer priming bias with VLMM --
    vlmm.fivep <- fitpar[[1]][["vlmm.fivep"]] #5'side
    vlmm.threep <- fitpar[[1]][["vlmm.threep"]] #3'side
    stopifnot(!is.null(vlmm.fivep))
    stopifnot(!is.null(vlmm.threep))
    
    #Calculate log(bias) for each fragment based on the VLMM
    fragtypes.temp <- addVLMMBias(fragtypes.temp, vlmm.fivep, vlmm.threep)
    
    # -- fit models --
    res[[1]] <- list()
    
    # remove first and last bp for predicting coverage along transcript
    not.first.or.last.bp <- !(fragtypes.temp$start == 1 | fragtypes.temp$end == l)
    fragtypes.temp <- fragtypes.temp[not.first.or.last.bp,]
    
    #Prepare the range/width of fragment read
    ir <- IRanges(fragtypes.temp$start, fragtypes.temp$end)
    
    #Output result
    res[[1]]$l <- l #Transcript length
    res[[1]]$frag.cov <- frag.cov #raw coverage
    res[[1]]$pred.cov <- list()
    
    gstart <- GenomicRanges::start(range(gene))
    gend <- GenomicRanges::end(range(gene))
    trx.length <- sum(width(ebt2[[test_ELAVL1_RXID]]))
    
    #for (modeltype in names(models)) {
    # message("predicting model type: ",modeltype)
    log.lambda <- getLogLambda(fragtypes.temp, models, modeltype, fitpar, bamfile)
    pred0 <- exp(log.lambda)
    pred <- pred0/mean(pred0)*mean(fragtypes.temp$count)
    #res[[1]][["pred.cov"]][[modeltype]] <- coverage(ir, weight=pred)
    
    fragtypes.temp.test <- fragtypes.temp
    fragtypes.temp.test$pred <- pred
    fragtypes.temp.test$ir <- ir
    fragtypes.temp.test2 <- fragtypes.temp.test[fragtypes.temp.test$count > 0,]
    #print(length(fragtypes.temp.test$count))
    #print(length(fragtypes.temp.test2$count))
    
    #pred.mapped <- fragtypes.temp.test2$pred*((length(fragtypes.temp.test$count)-length(fragtypes.temp.test2$count))/length(fragtypes.temp.test2$count))
    #ir.mapped <- fragtypes.temp.test2$ir
    pred.mapped <- c(0, fragtypes.temp.test2$pred*((length(fragtypes.temp.test$count)-length(fragtypes.temp.test2$count))/length(fragtypes.temp.test2$count)), 0)
    ir.mapped <- c(IRanges(1, 2), fragtypes.temp.test2$ir, IRanges(trx.length-1, trx.length))
    
    #print(pred.mapped)
    #print(ir.mapped)
    
    res[[1]][["pred.cov"]][["GC"]] <- coverage(ir.mapped, weight=pred.mapped)
    #}
    
    return(res)
}


#modeltype <- names(models)

#fragtypes <- fragtypes.temp
#models <- models
#modeltype <- modeltype[2]
#fitpar <- fitpar
#bamfile <- bamfile

getLogLambda <- function(fragtypes, models, modeltype, fitpar, bamfile) {
    #Formula
    f <- models[[modeltype]]$formula
    
    #Initialize offset
    offset <- numeric(nrow(fragtypes))
    
    ## -- fragment bias --
    if ("fraglen" %in% models[[modeltype]]$offset) {
        offset <- offset + fragtypes$logdfraglen
    }
    
    ## -- random hexamer priming bias with VLMM --
    if ("vlmm" %in% models[[modeltype]]$offset) {
        offset <- offset + fragtypes$fivep.bias + fragtypes$threep.bias
    }
    
    #Parameters
    if (!is.null(f)) {
        gc.knots <- seq(from = .4, to = .6, length = 3)
        gc.bk <- c(0, 1)
        relpos.knots <- seq(from = .25, to = .75, length = 3)
        relpos.bk <- c(0, 1)
        
        #Design matrix
        mm.big <- model.matrix(formula(f), data = fragtypes)
        
        #Parameters
        beta <- fitpar[[1]]$coef[[modeltype]]
        
        #Checking
        stopifnot(any(colnames(mm.big) %in% names(beta))) #Parameters exist or not.
        if (all(is.na(beta))) stop("all coefs are NA") #Parameters exist or not.
        
        #Replace NA coefs with 0: these were not observed in the training data
        beta[is.na(beta)] <- 0
        
        #This gets rid of the gene1, gene2... and Intercept terms
        beta <- beta[match(colnames(mm.big), names(beta))]
        
        #Add offset
        log.lambda <- as.numeric(mm.big %*% beta) + offset
    } else {
        log.lambda <- offset
    }
    #Checking
    if (!all(is.finite(log.lambda))) stop("log.lambda is not finite")
    
    return(log.lambda)
    
}

