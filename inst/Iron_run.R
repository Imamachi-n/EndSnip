############################################################
### Iron: Reducing fragment sequence bias in RNA Sequencing coverage.
### Created by Naoto Imamachi on 2016-02-24.
### Updated and maintained by Naoto Imamachi since Feb 2016.
############################################################

### This script was modified from alpine <https://github.com/mikelove/alpine>
###
### This program is free software; you can redistribute it and/or
### modify it under the terms of the GNU General Public License
### as published by the Free Software Foundation; either version 2
### of the License, or (at your option) any later version.
###
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
### GNU General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program; if not, see <http://www.gnu.org/licenses/>.
###
### 2016/02/24 Changed XXX....

#Required libraries
library(GenomicAlignments)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg19)
library(splines) # for GLM
library(speedglm) # for GLM

#Input filepath
#bamfile <- "/home/akimitsu/Documents/data/CFIm25_study/RNA-seq/COAD-Tumor-TCGA-A6-2675-01A-02R-1723-07/tophat_out/accepted_hits_chr19.bam"
#gtffile <- "/home/akimitsu/Documents/database/annotation_file/Refseq_gene_hg19_June_02_2014.gtf"
bamfile <- "C:/Users/Naoto/Documents/Visual Studio 2015/Projects/EndClip/EndClip/inst/data/accepted_hits_chr19.bam"
gtffile <- "C:/Users/Naoto/Documents/Visual Studio 2015/Projects/EndClip/EndClip/inst/data/Refseq_gene_hg19_June_02_2014_chr19.gtf"
singleUTRList <- "C:/Users/Naoto/Documents/Visual Studio 2015/Projects/EndClip/EndClip/DaPars_Test_data/EndClip_TCGA_Test_data_result_temp_extracted_chr19.txt"

#Read BAM file
ga <- readGAlignments(bamfile)
range(ranges(ga))

#Read gtf file
#circ_seqs: A character vector to list out which chromosomes should be marked as circular.
txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs = character())

#Prepare single-3UTR lists
su <- read.table(singleUTRList)
trxids <- as.vector(su[,1])

#Extract reference trxdb
txdf <- AnnotationDbi::select(txdb, 
                              columns = c("TXNAME","TXID", "GENEID"), 
                              keys = trxids, 
                              keytype = "TXNAME") #WARNING: which select function ?
txdf.txid <- txdf$TXID
txdf.txid <- txdf$TXID[1:20] #TEST: selected 10 transcripts with single-3UTR

#Extract exon information from gtf file
ebt <- exonsBy(txdb, "tx")

#Extract single-3UTR exon information
ebt <- ebt[txdf.txid]

#Prepare genome infor etc...
seqlevelsStyle(Hsapiens) <- "UCSC"
genenames <- names(ebt)
genenames <- genenames[1:20] #TEST: selected 10 transcripts with single-3UTR
names(genenames) <- genenames

#Load scripts
source("C:/Users/Naoto/Documents/Visual Studio 2015/Projects/EndClip/EndClip/inst/Iron-convertion_func.R")
source("C:/Users/Naoto/Documents/Visual Studio 2015/Projects/EndClip/EndClip/inst/Iron-data_prep_func.R")
source("C:/Users/Naoto/Documents/Visual Studio 2015/Projects/EndClip/EndClip/inst/Iron-VLMM.R")
source("C:/Users/Naoto/Documents/Visual Studio 2015/Projects/EndClip/EndClip/inst/Iron-GLM.R")

#Prepare two demensional matrix stored the information about the aligned paired-end fragments
#fragtypes <- lapply(genenames, function(gene) {
#    buildFragtypesFromExons(ebt[[gene]], genome = Hsapiens,
#                            readlength = 48, minsize = 100, maxsize = 300)
#})

# Make index file for BAM file (If bam.bai file does not exist, ...)
indexBam(bamfile)

# Prepare RNA-seq fragment sequence bias model
models <- list("GC" = list(formula = "count~ns(gc, knots = gc.knots, Boundary.knots = gc.bk) + gene",
                           offset = c("fraglen", "vlmm")))

## -- GLM fitting --
fitpar <- fitModelOverGenes(ebt, bamfile, fragtypes, Hsapiens, models,
                            readlength = 48, zerotopos = 2, speedglm = TRUE,
                            minsize = 100, maxsize = 300)

fitpar <- list(fitpar)
#names(fitpar) <- names(bamfiles)[1]

## -- Estimate unbiased sequence coverage --
ebt2 <- exonsBy(txdb, by = "tx")
trxnames2 <- names(ebt2)

#Extract TXNAME/GENEID/TXID data
geneids2 <- AnnotationDbi::keys(txdb, "TXNAME")
txdf2 <- AnnotationDbi::select(txdb, columns = c("TXNAME", "TXID", "GENEID"), keys = geneids2, keytype = "TXNAME")

#Fitting model
models <- list("null" = list(formula = NULL, offset = NULL),
               "GC" = list(formula = "count ~ ns(gc, knots = gc.knots, Boundary.knots = gc.bk) + 0",
               offset = c("fraglen", "vlmm")))

#TEST: 
test_ELAVL1_RXID <- txdf2[txdf2$GENEID == "ELAVL1",]$TXID
ebt2[[test_ELAVL1_RXID]]

test_res <- predictOneGene(ebt2[[test_ELAVL1_RXID]], bamfile, fitpar, genome=Hsapiens,
                           models, readlength = 48, minsize = 100, maxsize = 300)

gene <- ebt2[[test_ELAVL1_RXID]]
genome <- Hsapiens
readlength <- 48
minsize <- 100
maxsize <- 300

#Checking
stopifnot(is(gene, "GRanges"))
stopifnot(all(sapply(models, function(x) names(x) %in% c("formula","offset"))))
stopifnot(!is.null(fitpar))

bamfile <- list.files(bamfile)

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

## -- fragment bias --
fraglen.density <- fitpar[[1]][["fraglen.density"]]
stopifnot(!is.null(fraglen.density)) #Checking

#Calculate the density of fragment length from mapped reads
fragtypes.temp$logdfraglen <- log(matchToDensity(fragtypes.temp$fraglen, fraglen.density))

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

for (modeltype in names(models)) {
    # message("predicting model type: ",modeltype)
    log.lambda <- getLogLambda(fragtypes.temp, models, modeltype, fitpar, bamfile)
    pred0 <- exp(log.lambda)
    pred <- pred0/mean(pred0)*mean(fragtypes.temp$count)
    res[[1]][["pred.cov"]][[modeltype]] <- coverage(ir, weight=pred)
}

plotCov <- function(res, m="GC", cond, xlab="", ylab="", log=FALSE, lwd=3, ...) {
    for (i in seq_along(res)) {
        transform <- if (log) {
            function(x) log10(x + 1)
        } else {
            I
        }
        #ymax <- transform(getYmax(res))
        plot(transform(as.numeric(res[[i]]$frag.cov)), type="l", # ylim=c(0,ymax),
             xlab=xlab, ylab=ylab, col=cond[i], lwd=lwd, ...)
        lines(transform(as.numeric(res[[i]][["pred.cov"]][[m]])), col=rgb(0,0,0,.7),lwd=lwd)
    }
}

plotCov(res, cond="red")






