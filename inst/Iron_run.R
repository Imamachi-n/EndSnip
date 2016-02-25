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

#Extract exon information from gtf file
ebt <- exonsBy(txdb, "tx")

#Extract single-3UTR exon information
ebt <- ebt[txdf.txid]

#Prepare genome infor etc...
seqlevelsStyle(Hsapiens) <- "UCSC"
genenames <- names(ebt)
names(genenames) <- genenames

#TEST: buildFragtypesFromExons
exons <- ebt[[genenames[1]]]
genome <- Hsapiens
readlength <- 48
minsize <- 100
maxsize <- 300
npre <- 8
npost <- 12
gc <- TRUE
gc.str <- TRUE
vlmm <- TRUE

is(exons, "GRanges")
is(genome, "BSgenome")
is.numeric(minsize) & is.numeric(maxsize) & is.numeric(readlength)
sum(width(exons)) > maxsize
all(c("exon_rank", "exon_id") %in% names(mcols(exons)))





