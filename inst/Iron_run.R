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

print(paste0("[", Sys.time(), "] ", "Beginning Iron run (v0.1.0)"))
print("--------------------------------------------------")

#Required libraries
library(GenomicAlignments)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg19)
library(splines) # for GLM
library(speedglm) # for GLM

#Load scripts
source("C:/Users/Naoto/Documents/Visual Studio 2015/Projects/EndClip/EndClip/inst/Iron-convertion_func.R")
source("C:/Users/Naoto/Documents/Visual Studio 2015/Projects/EndClip/EndClip/inst/Iron-data_prep_func.R")
source("C:/Users/Naoto/Documents/Visual Studio 2015/Projects/EndClip/EndClip/inst/Iron-VLMM.R")
source("C:/Users/Naoto/Documents/Visual Studio 2015/Projects/EndClip/EndClip/inst/Iron-GLM.R")
source("C:/Users/Naoto/Documents/Visual Studio 2015/Projects/EndClip/EndClip/inst/Iron-predict_coverage.R")

#Input filepath
#bamfile <- "/home/akimitsu/Documents/data/CFIm25_study/RNA-seq/COAD-Tumor-TCGA-A6-2675-01A-02R-1723-07/tophat_out/accepted_hits_chr19.bam"
#gtffile <- "/home/akimitsu/Documents/database/annotation_file/Refseq_gene_hg19_June_02_2014.gtf"
bamfile <- "C:/Users/Naoto/Documents/Visual Studio 2015/Projects/EndClip/EndClip/inst/data/accepted_hits_chr10_NoCTRL_SE36.bam"
#bamfile <- "C:/Users/Naoto/Documents/Visual Studio 2015/Projects/EndClip/EndClip/inst/data/accepted_hits_chr19.bam"
gtffile <- "C:/Users/Naoto/Documents/Visual Studio 2015/Projects/EndClip/EndClip/inst/data/Refseq_gene_hg19_June_02_2014_chr10.gtf"
singleUTRList <- "C:/Users/Naoto/Documents/Visual Studio 2015/Projects/EndClip/EndClip/DaPars_Test_data/EndClip_TCGA_Test_data_result_temp_extracted_chr10.txt"

#Bam file information
read.type <- "SE" #PE
read.length <- 36

#Read BAM file
ga.mapped <- readGAlignments(bamfile)
range(ranges(ga.mapped))

# Make index file for BAM file (If bam.bai file does not exist, ...)
indexBam(bamfile)

#Read gtf file
#circ_seqs: A character vector to list out which chromosomes should be marked as circular.
txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs = character()) #Core data

#Extract exon information from gtf file
ebt <- exonsBy(txdb, "tx") #Core data

#Prepare single-3UTR lists
su <- read.table(singleUTRList)
trxids <- as.vector(su[,1])

#Extract reference trxdb
txdf <- AnnotationDbi::select(txdb, 
                              columns = c("TXNAME","TXID", "GENEID"), 
                              keys = trxids, 
                              keytype = "TXNAME") #WARNING: which select function ?
txdf.txid <- txdf$TXID
txdf.txid <- txdf$TXID[1:70] #TEST: selected 10 transcripts with single-3UTR

#Re-extract reference trxdb using gene symbol
txdf.geneid <- unique(txdf$GENEID)
txdf.re <- AnnotationDbi::select(txdb,
                                 columns = c("TXNAME", "TXID", "GENEID"),
                                 keys = txdf.geneid,
                                 keytype = "GENEID")
txdf.re.txid <- txdf.re$TXID #All isoform with single-UTR
txdf.re.geneid <- unique(txdf.re$GENEID) #All gene names with single-UTR

#Define representative isoforms
trxid.rep.list <- c()
for (geneid in txdf.re.geneid) {
    #geneid <- txdf.re.geneid[4] #TEST:
    #Extract all isoforms for each gene
    curr.geneid.trxid <- txdf.re[txdf.re$GENEID == geneid,]$TXID
    
    if (length(curr.geneid.trxid) == 1) {
        trxid.rep.list <- append(trxid.rep.list, curr.geneid.trxid)
        next
    }
    
    curr.geneid.test <- list()
    for (trxid in curr.geneid.trxid) {
        #trxid <- curr.geneid.trxid[1] #TEST:
        
        #The range of test_gene (GRanges object)
        gene <- ebt[[trxid]]
        generange <- range(gene)
        strand(generange) <- "*" #Not necessary
        
        #Mapped reads on each transcript (st -> ed (including exon/intron))
        suppressWarnings({
            ga <- readGAlignAlpine(bamfile, generange, read.type) #generange should be GRanges (Not GRangeList!!)
        })
        
        #Single-end TEST:
        #ga <- findOverlaps(ga.mapped, ebt[trxid], ignore.strand = TRUE)
        
        #nhitPerQuery <- function(x) tabulate(queryHits(x), nbins = queryLength(x))
        #ga.ntx <- nhitPerQuery(ga)
        #mcols(ga.mapped)$ntx <- ga.ntx
        #table(ga.ntx)
        
        #ga.grl <- grglist(ga.mapped, order.as.in.query = TRUE) #GAlignments -> GRangesList
        #ga.grlf <- flipQuery(ga.grl) #Flipped
        
        #ga <- readGAlignments(bamfile, param = ScanBamParam(which = generange, flag = scanBamFlag(isSecondaryAlignment = FALSE)))
        #fco <- findCompatibleOverlaps(ga, GRangesList(gene))
        
        ga <- keepSeqlevels(ga, seqlevels(generange))
        
        #Compatible reads on each transcript (st -> ed (including only exon!!))
        fco <- findCompatibleOverlaps(ga, GRangesList(gene))
        
        #Transcript position of compatible reads
        reads <- gaToReadsOnTx(ga, GRangesList(gene), fco, read.type, read.length)
        #grl <- GRangesList(gene)
        trx.length <- sum(width(gene))
        trx.exp <- length(reads[[1]])/trx.length*1000
        
        trx.cov <- c(0, as.vector(coverage(reads[[1]])), 0)
        
        #Result
        curr.geneid.test$trxid <- append(curr.geneid.test$trxid, trxid)
        curr.geneid.test$trx.exp <- append(curr.geneid.test$trx.exp, trx.exp)
        
        #curr.geneid.test$test <- c("test")
        #curr.geneid.test$test <- append(curr.geneid.test$test,"test")
    }
    max.index <- which.max(curr.geneid.test$trx.exp)
    trxid.rep <- curr.geneid.test$trxid[max.index]
    
    trxid.rep.list <- append(trxid.rep.list, trxid.rep)
}

txdf.rep <- AnnotationDbi::select(txdb,
                                  columns = c("TXNAME", "TXID", "GENEID"),
                                  keys = as.character(trxid.rep.list),
                                  keytype = "TXID")

ebt.rep <- ebt[as.character(trxid.rep.list)]

#Extract single-3UTR exon information
#ebt <- ebt[txdf.txid]

#Prepare genome infor etc...
seqlevelsStyle(Hsapiens) <- "UCSC"
genenames <- names(ebt)
#genenames <- genenames[1:20] #TEST: selected 10 transcripts with single-3UTR
names(genenames) <- genenames

#Prepare two demensional matrix stored the information about the aligned paired-end fragments
#fragtypes <- lapply(genenames, function(gene) {
#    buildFragtypesFromExons(ebt[[gene]], genome = Hsapiens,
#                            readlength = 48, minsize = 100, maxsize = 300)
#})

# Prepare RNA-seq fragment sequence bias model
models <- list("GC" = list(formula = "count~ns(gc, knots = gc.knots, Boundary.knots = gc.bk) + gene",
                           offset = c("fraglen", "vlmm")))

#SE_version
models <- list("GC" = list(formula = "count~ns(gc, knots = gc.knots, Boundary.knots = gc.bk) + gene",
                              offset = c("vlmm")))

#models <- list("GC" = list(formula = "count~ns(relpos, knots = relpos.knots, Boundary.knots = relpos.bk) + gene",
#                           offset = c("fraglen", "vlmm")))


#models <- list("GC" = list(formula = "count~ns(gc, knots = gc.knots, Boundary.knots = gc.bk) + gene",
#                           offset = c("vlmm")))

#ebt.last <- GRangesList(lapply(ebt, function(x){
#    x[length(x)]
#}))

## -- GLM fitting --
fitpar <- fitModelOverGenes(ebt.rep[1:70], bamfile, Hsapiens, models, read.type,
                            readlength = read.length, zerotopos = 2, speedglm = TRUE,
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
#models <- list("null" = list(formula = NULL, offset = NULL),
#               "GC" = list(formula = "count ~ ns(gc, knots = gc.knots, Boundary.knots = gc.bk) + 0",
#               offset = c("fraglen", "vlmm")))

models <- list("null" = list(formula = NULL, offset = NULL),
               "GC" = list(formula = "count ~ ns(gc, knots = gc.knots, Boundary.knots = gc.bk) + 0",
               offset = c("vlmm")))

#models <- list("null" = list(formula = NULL, offset = NULL),
#               "GC" = list(formula = "count ~ ns(relpos, knots = relpos.knots, Boundary.knots = relpos.bk) + 0",
#                           offset = c("fraglen", "vlmm")))

#models <- list("null" = list(formula = NULL, offset = NULL),
#               "GC" = list(formula = "count ~ ns(gc, knots = gc.knots, Boundary.knots = gc.bk) + 0",
#                           offset = c("vlmm")))


#TEST: ELAVL1
test_ELAVL1_RXID <- txdf2[txdf2$GENEID == "ELAVL1",]$TXID
map2 <- mapTxToGenome(ebt2[[test_ELAVL1_RXID]])

res <- predictOneGene(ebt2[[test_ELAVL1_RXID]], bamfile, fitpar, genome=Hsapiens,
                      models, readType, readlength = 48, minsize = 100, maxsize = 300)

#TEST: ELAVL1
map2$cov <- as.vector(res[[1]]$pred.cov$GC)
chrom_number <- seqlevels(ebt2[[test_ELAVL1_RXID]])
bedgraph <- data.frame(chrom=rep(chrom_number,sum(width(ebt2[[test_ELAVL1_RXID]]))), st=map2$genome-1, ed=map2$genome, cov=map2$cov)
write.table(bedgraph, file="chr19_ELAVL1.bg", quote=F, sep="\t", row.names=F, col.names=F)

#Checking
plotCov <- function(res, m="GC", cond, xlab="", ylab="", log=FALSE, lwd=3, ...) {
    for (i in seq_along(res)) {
        transform <- if (log) {
            function(x) log10(x + 1)
        } else {
            I
        }
        #ymax <- transform(getYmax(res))
        plot(transform(as.numeric(res[[i]]$frag.cov)), type="l",  ylim=c(0,max(as.vector(res[[1]]$pred.cov$GC), as.vector(res[[1]]$frag.cov))),
             xlab=xlab, ylab=ylab, col=cond[i], lwd=lwd, ...)
        lines(transform(as.numeric(res[[i]][["pred.cov"]][[m]])), col=rgb(0,0,0,.7),lwd=lwd, ylim=c(0,max(as.vector(res[[1]]$pred.cov$GC))))
    }
}

plotCov(res, cond="red")






