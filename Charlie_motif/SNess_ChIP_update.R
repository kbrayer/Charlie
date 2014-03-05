##rGADEM from a BED file

##starting from BED files

bed <- "bedfile.bed"

##libraries
library(ChIPpeakAnno)
data(TSS.human.NCBI36)

##Import files
a <- read.table(bed, header=FALSE, skip=2, sep="\t", row.names=NULL, 
                col.names=c("chr", "start", "end", "name", "score", "strand", "NA", "ID?"))
matrixa <- as.matrix(a)
length(a[,1])

a.bed = data.frame(matrixa)

a.bed[,6] <- 1 ##Set strand to 1 ??

## convert to Range data ##

a.rangedData = BED2RangedData(a.bed)

##get annotations##
a.annotatedPeak = annotatePeakInBatch(a.rangedData, AnnotationData = TSS.human.NCBI36)

##get summary info##
pie(table(as.data.frame(a.annotatedPeak)$insideFeature))
saved.pie <- recordPlot()

##problem starting here
filename <- paste(name1, "_pie.jpg", sep="", collapse="")
jpeg(filename, height = 700, width = 700)
replayPlot(saved.pie)
dev.off()
