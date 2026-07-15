> tabla <- read.xls("ChiPanneW.xls")
> colnames(tabla)
 [1] "Chr"                "Start"              "Stop"              
 [4] "length"             "id"                 "pval"              
 [7] "X.10..log10.Pval.." "peak"               "X"                 
[10] "X.1"                "X.2"                "X.3"               
> tabla <- tabla[,c(1:8)]
> colnames(tabla)
[1] "Chr"                "Start"              "Stop"              
[4] "length"             "id"                 "pval"              
[7] "X.10..log10.Pval.." "peak"              
> head(tabla)
   Chr   Start    Stop length id         pval X.10..log10.Pval..     peak
1 chr1 4774743 4776846   2103  5 7.182253e-04           31.43739 Increase
2 chr1 4796955 4798978   2023  6 2.990000e-11          105.24688 Increase
3 chr1 4846900 4849212   2312  7 6.599923e-03           21.80461 Increase
4 chr1 5072328 5074300   1972 11 9.320000e-06           50.30609 Increase
5 chr1 6203703 6206190   2487 13 2.340000e-06           56.30310 Increase
6 chr1 6296014 6297441   1427 14 1.452488e-03           28.37887 Increase

> library(GenomicRanges)
> require(GenomicRanges)
> head(tabla)
   Chr   Start    Stop length id         pval X.10..log10.Pval..     peak
1 chr1 4774743 4776846   2103  5 7.182253e-04           31.43739 Increase
2 chr1 4796955 4798978   2023  6 2.990000e-11          105.24688 Increase
3 chr1 4846900 4849212   2312  7 6.599923e-03           21.80461 Increase
4 chr1 5072328 5074300   1972 11 9.320000e-06           50.30609 Increase
5 chr1 6203703 6206190   2487 13 2.340000e-06           56.30310 Increase
6 chr1 6296014 6297441   1427 14 1.452488e-03           28.37887 Increase
> Chip <- GRanges(tabla$Chr,IRanges(tabla$Start,tabla$Stop),strand='*')
> names(Chip) <- paste("AWpeakZic2KD",1:nrow(tabla),sep="_")
> colnames(tabla)
[1] "Chr"                "Start"              "Stop"              
[4] "length"             "id"                 "pval"              
[7] "X.10..log10.Pval.." "peak"              
> score(Chip) <- tabla[,"X.10..log10.Pval.."]
> Chip  ## Show the first and last lines of the GRanges object
	
GRanges object with 33143 ranges and 1 metadata column:
                     seqnames                 ranges strand |            score
                        <Rle>              <IRanges>  <Rle> |        <numeric>
      AWpeakZic2KD_1     chr1     [4774743, 4776846]      * | 31.4373932824388
      AWpeakZic2KD_2     chr1     [4796955, 4798978]      * | 105.246880241934
      AWpeakZic2KD_3     chr1     [4846900, 4849212]      * |  21.804611575056
      AWpeakZic2KD_4     chr1     [5072328, 5074300]      * | 50.3060859595968
      AWpeakZic2KD_5     chr1     [6203703, 6206190]      * | 56.3030984231021
                 ...      ...                    ...    ... .              ...
  AWpeakZic2KD_33139     chrX [164868966, 164870101]      * | 72.0609567575279
  AWpeakZic2KD_33140     chrX [165032917, 165033890]      * | 20.3584465308559
  AWpeakZic2KD_33141     chrX [165233134, 165234213]      * | 18.2846243352497
  AWpeakZic2KD_33142     chrX [166276021, 166276848]      * | 29.1671536526436
  AWpeakZic2KD_33143     chrX [166435012, 166436274]      * | 66.4554094470433
  -------
  seqinfo: 20 sequences from an unspecified genome; no seqlengths

################################################################################
####Basic number and statistics of the Chip peaks

##3.1 What is the number of peaks?
> length(Chip)
[1] 33143


> Chip.size <- width(Chip)  ## What is the mean, median, and max size of the peaks?
				## The function width() returns the sizes of all ranges in a GRanges object as vector. 
> summary(Chip.size)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    488    1508    2139    2938    3302   63292 

##Remove all peaks that are larger than 2kb
#A GRanges object can be indexed using the [] operator similar than vectors. 

Chip1 <- Chip[width(Chip)<=2000]

################################################################################

###What is the distribution of peak sizes?

#Plot the distribution of sizes using the hist() function.

pdf("distributionofsizes_ChipAW")
hist(width(Chip1), xlab="Chip Zic2KOD (AW) peak size", col="gray")
dev.off()
################################################################################

###What is the distribution of  ChIP-seq peak p-values?

#We plot the histogram of −log10 transformed p-values. Remember that we stored them in the score column of the GRanges object. 
#We can access them using the function score(). 

mlPvals <- score(Chip1)

pdf("distributionofChIPseqpeak_pvalues")
hist(mlPvals, xlab="-log_10(p-value)", col="gray")
dev.off()

################################################################################
# Compare the peaks of ChipEH and chipAW
> tablaEH <- read.xls("ChIpEH.xls")
> ChipEH <- GRanges(tablaEH$Chr,IRanges(tablaEH$Start,tablaEH$Stop),strand='*')
names(ChipEH) <- paste("peakZic2OE",1:nrow(tablaEH),sep="_")
colnames(tablaEH)
[1] "Chr"                "Start"              "Stop"              
[4] "Width"              "id"                 "Pval"              
[7] "X.10..log10.Pval.."
> score(ChipEH) <- tablaEH[,"X.10..log10.Pval.."]
ChipEH

Chip1EH <- ChipEH[width(ChipEH)<=2000]

################################################################################
##  Make a barplot of the number of peaks

pdf("barplot_numberofpeaks")
bp <- barplot(c(length(Chip1EH), length(Chip1)), names=c("EH", "AW"))
# add actual values as text lables to the plot
text(bp, c(length(Chip1EH), length(Chip1)), labels=c(length(Chip1EH), length(Chip1)), pos=1)
dev.off()

################################################################################

# How many EH peaks overlap ChipAW peaks?
#We can calculate the overlaps between two sets of regions with the findOverlaps function from the GenomicRanges package. 
#To get the subset of EH peaks that overlap with the ChipAW peaks we can use the ‘subsetByOverlaps’ function.

# compute overlaps
ovlHits <- findOverlaps(Chip1EH, Chip1)

# show the resulting Hits object
ovlHits

#Hits object with 0 hits and 0 metadata columns:
#   queryHits subjectHits
#   <integer>   <integer>
#  -------
#  queryLength: 4435 / subjectLength: 15093

# get subsets of binding sites
ovl <- subsetByOverlaps(Chip1EH, Chip1)

# number of overlaps
length(ovl)

# as percent
length(ovl) / length(Chip1EH) * 100

#  No hay picos en comun
################################################################################
# Make a Venn-diagram showing the overlap of of ER and FOXA1 peaks.
# We first take the subset of peaks that are unique to EH and AW using the function setdiff.
# get subsets of binding sites
EH.uniq <- setdiff(Chip1EH, Chip1)
AW.uniq <- setdiff(Chip1, Chip1EH)

# Then we use the package Vennerable to build a Venn object with the number of peaks in each subset. This object can than be passed to the ‘plot’ function.
# plot overlap of binding sites as Venn diagram
require(Vennerable)

# build objects with the numbers of sites in the subsets
venn <- Venn(SetNames=c("Chip1EH", "Chip1"), 
    Weight=c(
        '10'=length(EH.uniq), 
        '11'=length(ovl), 
        '01'=length(AW.uniq)
    )
)

# plot Venn Diagram
pdf("venny_Chipseqpeaks.pdf")
plot(venn)
dev.off()
################################################################################

## Functional annotation of ChIP-seq peaks

#To understand the function of a transcription factor we want to know to which genes and genomic features it binds.

#We will use a TxDb objects. Such an object is an R interface to prefabricated databases contained by specific annotation packages. 
#The package TxDb.Mmusculus.UCSC.mm10.ensGene includes all human genes and transcripts from UCSC with coordinates for the mm10 genome assembly.

source("https://bioconductor.org/biocLite.R")
biocLite("TxDb.Mmusculus.UCSC.mm10.ensGene")
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene

Now we want to assign each peak to a chromosomal region by using the feature data base.
require(ChIPpeakAnno)   # laod package to annotate peaks with genes

biocLite("ChIPpeakAnno")
library(dplyr)
library(ChIPpeakAnno)
# calucluate the overlap with features
EH.features <- assignChromosomeRegion(Chip1EH, TxDb=txdb, nucleotideLevel=FALSE)
# show the results
EH.features




################################################################################







