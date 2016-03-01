predictOneGene <- function(gene, bamfile, fitpar, genome=Hsapiens,
                           models, readlength, minsize, maxsize) {
    #Checking
    stopifnot(is(gene, "GRanges"))
    stopifnot(all(sapply(models, function(x) names(x) %in% c("formula","offset"))))
    stopifnot(!is.null(fitpar))

    fragtypes <- buildFragtypesFromExons(gene, genome, readlength=readlength,
                                         minsize=minsize, maxsize=maxsize)
    res <- list()
    bamname <- bamfile
    # add counts
    generange <- range(gene)
    strand(generange) <- "*" # not necessary
    suppressWarnings({
        ga <- readGAlignAlpine(bamfile, generange)
    })
    if (length(ga) == 0) {
        res[[bamname]] <- as.list(rep(NA,length(models)))
        names(res[[bamname]]) <- names(models)
        next
    }
    ga <- keepSeqlevels(ga, as.character(seqnames(gene)[1]))
    fco <- findCompatibleOverlaps(ga, GRangesList(gene))
    # message("-- ",round(length(fco)/length(ga),2)," compatible overlaps")
    reads <- gaToReadsOnTx(ga, GRangesList(gene), fco)
    
    # save fragment coverage for later
    l <- sum(width(gene))
    frag.cov <- coverage(reads[[1]][start(reads[[1]]) != 1 & end(reads[[1]]) != l])
    
    fragtypes.temp <- matchReadsToFraglist(reads, list(fragtypes))[[1]]
    ## -- fragment bias --
    fraglen.density <- fitpar[[bamname]][["fraglen.density"]]
    stopifnot(!is.null(fraglen.density))
    fragtypes.temp$logdfraglen <- log(matchToDensity(fragtypes.temp$fraglen, fraglen.density))
    ## -- random hexamer priming bias with VLMM --
    vlmm.fivep <- fitpar[[bamname]][["vlmm.fivep"]]
    vlmm.threep <- fitpar[[bamname]][["vlmm.threep"]]
    stopifnot(!is.null(vlmm.fivep))
    stopifnot(!is.null(vlmm.threep))
    fragtypes.temp <- addVLMMBias(fragtypes.temp, vlmm.fivep, vlmm.threep)
    
    # -- fit models --
    res[[bamname]] <- list()
    
    # remove first and last bp for predicting coverage along transcript
    not.first.or.last.bp <- !(fragtypes.temp$start == 1 | fragtypes.temp$end == l)
    fragtypes.temp <- fragtypes.temp[not.first.or.last.bp,]
    
    ir <- IRanges(fragtypes.temp$start, fragtypes.temp$end)
    res[[bamname]]$l <- l
    res[[bamname]]$frag.cov <- frag.cov
    res[[bamname]]$pred.cov <- list()
    for (modeltype in names(models)) {
        # message("predicting model type: ",modeltype)
        log.lambda <- getLogLambda(fragtypes.temp, models, modeltype, fitpar, bamname)
        pred0 <- exp(log.lambda)
        pred <- pred0/mean(pred0)*mean(fragtypes.temp$count)
        res[[bamname]][["pred.cov"]][[modeltype]] <- coverage(ir, weight=pred)
    }
res    
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

