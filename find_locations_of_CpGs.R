

### https://www.r-bloggers.com/using-genomation-to-analyze-methylation-profiles-from-roadmap-epigenomics-and-encode/


Using genomation to analyze methylation profiles from Roadmap epigenomics and ENCODE

The genomation package is a toolkit for annotation and visualization of various genomic data. The package is currently in developmental version of BioC. 
It allows to analyze high-throughput data, including bisulfite sequencing data. 
Here, we will visualize the distribution of CpG methylation around promoters and their locations within gene structures on human chromosome 3.

################################################################################################################################################################################################
Heatmap and plot of meta-profiles of CpG methylation around promoters

In this example we use data from Reduced Representation Bisulfite Sequencing (RRBS) and
Whole-genome Bisulfite Sequencing (WGBS) techniques and H1 and IMR90 cell types
derived from the ENCODE and the Roadmap Epigenomics Project databases.

We download the datasets and convert them to GRanges objects. Using rtracklayer and genomation functions. We also use a refseq bed file for annotation and extraction of promoter regions using readTranscriptFeatures function.

# download transcript data
refseq.path <- tempfile()
download.file("https://dl.dropboxusercontent.com/u/1373164/chr3.refseq.hg19.bed",destfile=refseq.path, method="curl")


# get RRBS from ENCODE
rrbs.IMR90.path <- tempfile()
download.file("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibMethylRrbs/wgEncodeHaibMethylRrbsImr90UwSitesRep1.bed.gz",
              ,destfile=rrbs.IMR90.path, method="curl")
rrbs.H1.path <- tempfile()
download.file("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibMethylRrbs/wgEncodeHaibMethylRrbsH1hescHaibSitesRep1.bed.gz",
              ,destfile=rrbs.H1.path, method="curl")

# read the RRBS data to GRanges
library(genomation)

rrbs.IMR90<-readGeneric(rrbs.IMR90.path, chr = 1, start = 2, end = 3, strand = 6,
            meta.cols = list(score=11), keep.all.metadata = FALSE, zero.based = T)
rrbs.H1<-readGeneric(rrbs.H1.path, chr = 1, start = 2, end = 3, strand = 6,
                        meta.cols = list(score=11), keep.all.metadata = FALSE, zero.based = T)



# mapping to [0, 1] range
rrbs.IMR90$score<-rrbs.IMR90$score/100
rrbs.H1$score<-rrbs.H1$score/100


# download WGBS data from Roadmap Epigenomics, convert to GRanges
library(rtracklayer)
wgbs.H1.path <- "http://egg2.wustl.edu/roadmap/data/byDataType/dnamethylation/WGBS/FractionalMethylation_bigwig/E003_WGBS_FractionalMethylation.bigwig"
wgbs.IMR90.path <- "http://egg2.wustl.edu/roadmap/data/byDataType/dnamethylation/WGBS/FractionalMethylation_bigwig/E017_WGBS_FractionalMethylation.bigwig"
grange <- GRanges("chr3", IRanges(1, 2e8)) #chr3
wgbs.H1 <- import.bw(wgbs.H1.path, selection=BigWigSelection(grange))
wgbs.IMR90 <- import.bw(wgbs.IMR90.path, selection=BigWigSelection(grange))



# Then we read transcripts from the refseq bed file 
# within 4 bp distance from the longest transcript of each gene.

transcriptFeat=readTranscriptFeatures(refseq.path, up.flank = 2000, down.flank = 2000)

################################################################################################################################################################################################

########################################################
#
# Heatmap produced from the analysis
#
########################################################

# Plot between 2 condition of the interest
d1=read.csv("C1_Paj.csv",header=T, sep="\t")
d2=read.csv("C2_Paj.csv",header=T, sep="\t")

# Merge data
data=merge(d1,d2,by="X",incomparables=NA,all=TRUE)
rn=rownames(data)
unique(rn)
rn=data[,1]
colnames(data)=c("Gene","C1","C2")
data2=sapply(data,function(x) if (is.factor(x)) { as.numeric(as.character(x))}else{x})
rownames(data2)=data[,1]
data2=data2[,2:3]
write.table(data2,"./order.csv",sep="\t")
hm <- heatmap.2(data2, scale="col", Rowv=F, Colv=F, symkey=FALSE,
                margins=c(8,8), cexRow=0.7, cexCol=1.0, key=TRUE, keysize=1.5,
                trace="none",density.info=c("density"),tracecol="blue",col=redgreen(100), main="C1 vs
C2")

# Example of the selection
slist=data[1350:1370,c(2,3)]
slist2=sapply(slist,function(x) if (is.factor(x)) { as.numeric(as.character(x))}else{x})
rownames(slist2)=data[1350:1370,1]
hm <- heatmap.2(slist2, scale="col", Rowv=F, Colv=F, symkey=FALSE,
                margins=c(8,8), cexRow=0.7, cexCol=1.0, key=TRUE, keysize=1.5,
                trace="none",density.info=c("density"),tracecol="blue",col=redgreen(100),
                main="FoldChange C1 vs. C2")


################################################################################################################################################################################################

Since we have read the files now we can build base-pair resolution matrices of scores(methylation values) for each experiment. The returned list of matrices can be used to draw heatmaps or meta profiles of methylation ratio around promoters. 

# Once we have read the files
# now we can build base-pair resolution matrices of scores for each experiment.
# The returned list of matrices can be used to draw meta profiles or heatmaps of read coverage around promoters.


targets <- list(WGBS.H1=wgbs.H1, WGBS.IMR90=wgbs.IMR90, RRBS.H1=rrbs.H1, RRBS.IMR90=rrbs.IMR90) 
sml <- ScoreMatrixList(targets=targets, windows=transcriptFeat$promoters, bin.num=20, strand.aware=TRUE, weight.col="score", is.noCovNA = TRUE)
 
sml.sub =intersectScoreMatrixList(sml, reorder = FALSE) # because scoreMatrices have different dimensions, we need them to cover the same promoters for the heatmap
multiHeatMatrix(sml.sub, xcoords=c(-2000, 2000), matrix.main=names(targets), xlab="region around TSS")


################################################################################################################################################################################################

smlcolors <- rainbow(4)
plotMeta(sml.sub, line.col=smlcolors)
legend("bottomright", 
       names(targets),
       lty=c(1,1), 
       lwd=c(2.5,2.5),
       col=smlcolors)

################################################################################################################################################################################################


Distribution of covered CpGs across gene regions
genomation facilitates visualization of given locations of features aggregated by  exons, introns, promoters and TSSs. To find the distribution of covered CpGs within these gene structures, we will use transcript features we previously obtained. Here is the breakdown of the code

1    Count overlap statistics between our CpGs from WGBS and RRBS H1 cell type and gene structures
# Finding locations of CpGs annotated by genic parts 

# The genomation facilitate visualization of given locations of features annotated by 
# exons, introns, promoters and TSSs.
# To find the distribution of CpGs around these gene structures, we will 
# first read the transcript features from a file using the readTranscriptFeatures function.

# Using genomation you can easly count overlap statistics between our CpGs from WGBS and RRBS and gene structures,
ann.wgbs <- annotateWithGeneParts(wgbs.H1, transcriptFeatures ) 
ann.rrbs <- annotateWithGeneParts(rrbs.H1, transcriptFeatures ) 

# calculate percentage of CpGs overlapping with annotation
getTargetAnnotationStats(ann.wgbs, percentage=TRUE, precedence=TRUE)
getTargetAnnotationStats(ann.rrbs, percentage=TRUE, precedence=TRUE)


2    Calculate percentage of CpGs overlapping with annotation plot them in a form of pie charts

# ..and plot them in a form of pie charts
par(mfrow=c(1,2))
plotTargetAnnotation(ann.wgbs, main="WGBS H1")
plotTargetAnnotation(ann.rrbs, main="RRBS H1")

################################################################################################################################################################################################

http://www.kumc.edu/Documents/biostatistics/Methylation%20workshop.pdf




################################################################################################################################################################################################
