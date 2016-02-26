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
txdf <- select(txdb, columns = c("TXNAME", "TXID", "GENEID"), keys = trxids, keytype = "TXNAME")
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

#Prepare two demensional matrix stored the information about the aligned paired-end fragments
#fragtypes <- lapply(genenames, function(gene) {
#    buildFragtypesFromExons(ebt[[gene]], genome = Hsapiens,
#                            readlength = 48, minsize = 100, maxsize = 300)
#})

# Make index file for BAM file (If bam.bai file does not exist, ...)
indexBam(bamfile)

# Prepare RNA-seq fragment sequence bias model
models <- list("GC" = list(formaula = "count~ns(gc, knots = gc.knots, Boundary.knots = gc.bk) + gene",
                           offset = c("fraglen", "vlmm")))

#TEST: fitModelOverGenes
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
readGAlignAlpine <- function(bamfile, generange){
    readGAlignmentPairs(bamfile, param = ScanBamParam(which = generange, flag = alpineFlag()))
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

genes <- ebt
bamfile <- bamfile
fragtypes <- fragtypes
genome <- Hsapiens
models <- models
readlength <- 48
zerotopos <- 2
speedglm <- TRUE
minsize <- 100
maxsize <- 300

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
        ga <- readGAlignAlpine(bamfile, generange)
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
    fco <- findCompatibleOverlaps(ga, GRangesList(gene))
    
    #Transcript position of compatible reads
    reads <- gaToReadsOnTx(ga, GRangesList(gene), fco)
    
    #Prepare dummy data(All possible fragment patterns)
    fragtypes <- list(buildFragtypesFromExons(ebt[[gene.name]], genome = Hsapiens,
                                              readlength = 48, minsize = 100, maxsize = 300))
    
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

#Check the numeber of fragments(mapped reads + dummy reads)
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

# -- Random hexamer priming bias with VLMM --







