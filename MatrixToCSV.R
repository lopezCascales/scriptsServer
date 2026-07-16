library(Matrix)
readLines("matrix.mtx")

#[17451] "36612 378 1"                                     
#[17452] "36624 378 1"                                     
#[17453] "36626 378 5"                                     
#[17454] "36634 378 1"                                     
#[17455] "36641 378 30"                                    
#[17456] "36643 378 1"                                     
#[17457] "36654 378 15"                                    
#[17458] "36659 378 2"                                     
#[17459] "36663 378 3"                                     
#[17460] "36669 378 26"                                    
#[17461] "36670 378 1"                                     
#[17462] "36672 378 1"                                     
#[17463] "36678 378 2"                                     
#[17464] "36693 378 1"  

mat <- readMM("matrix.mtx")

str(mat)
#Formal class 'dgTMatrix' [package "Matrix"] with 6 slots
#  ..@ i       : int [1:17461] 524 3113 5217 5765 6456 6680 8990 15338 19838 20734 ...
#  ..@ j       : int [1:17461] 0 0 0 0 0 0 0 0 0 0 ...
#  ..@ Dim     : int [1:2] 36693 378
#  ..@ Dimnames:List of 2
#  .. ..$ : NULL
#  .. ..$ : NULL
#  ..@ x       : num [1:17461] 1 1 2 4 1 1 1 1 1 1 ...
#  ..@ factors : list()

barcodes <- read.csv("barcodes.tsv", header=F)
#colnames(barcodes) <- c("barcodes")
features <- read.csv("features.tsv", sep="\t", header=F)
#colnames(features) <- c("ensemble","GENEID","type")


head(barcodes)
head(features)
colnames(mat) <- barcodes$V1
rownames(mat) <- features$V2

colSums(mat)
#AAACAGGC AAAGCGGA AAAGGCTG AACACGCA AACATGGG AACCCAAC AACCGCTT AACCGGAA 
#     180      205      185      198      177      215      168      288 
#AACCTGCT AACGAGGT AACTCTGG AAGACAGC AAGCACAT AAGCGAGT AAGCTGCA AAGGAGCA 
#     317      453      207      189      183      239      208      209 
#AAGGTCAG AAGTCTCG AAGTGGCT AATAGGGC AATCATGC AATCGGCA AATGGTGG ACAAACCG 
#     148      145      157      197      133      149      155      104 
#ACAAGGCA ACACCGTG ACAGAAGC ACAGAGCT ACAGGCAT ACATGTGC ACCAAGGA ACCACGAT 

head(mat)
df <- as.data.frame(as.matrix(mat))
write.csv(df,"TranscriptCounts.csv")
