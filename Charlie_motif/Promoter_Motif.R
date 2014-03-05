#An attempt to look at TF binding sites in promoter regions of significant genes from RNA-Seq Data
#Charlie Brayer Feb 2014

#This script will take output from RNA seq data analysis, obtain promoter sequences,
#look for specific motifs and score promoters - somehow. Hopefully will eventually also
#scan promoters for unknown TF binding sites. 

###Get promoter sequence for differentially expressed genes.

library(MotifDb)
library(seqLogo)
library(motifStack)
library(Biostrings)
library(GenomicFeatures)
library(hgu133plus2.db)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)


query (query (MotifDb, "Hsapiens"),"MYB")
#query(MotifDb, "MA0100.2")
query(MotifDb, "MYB")
pfm.MYB.jaspar <- query (query (MotifDb, "mmusculus"),"MYB") [[3]]
seqLogo(pfm.MYB.jaspar)

# Test with ABCB4

ABCB4 <- "5244"
chromosomal.loc <- transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by = "gene")[ABCB4]
promoter.ABCB4 <- getPromoterSeq(chromosomal.loc, Hsapiens, upstream = 1000, 
                                downstream = 0)
pcm.MYB.jaspar <- round(100 * pfm.MYB.jaspar)
matchPWM(pcm.MYB.jaspar, unlist(promoter.ABCB4)[[1]], "90%")

##Test with several genes


SigGenes <- read.table("Significant_RNA_seq_genes.txt", header=TRUE, sep ='\t')


Entrezgenes <- SigGenes[,3]
length(Entrezgenes)
head(Entrezgenes)

egREFSEQ <- toTable(org.Hs.egREFSEQ)
m <- match(SigGenes$ID, egREFSEQ$accession)
SigGenes$EntrezGene <- egREFSEQ$gene_id[m]
GetProms <- as.character(SigGenes [,8])
#sort(GetProms, decreasing = TRUE, na.last=NA)


#txdb2 <- makeTranscriptDbFromUCSC(genome="hg19", tablename="knownGene",
                                  #transcript_ids=SigGenes)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
grl <- transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by="gene")
#grl <- transcriptsBy(txdb, by="gene") 
fail <- setdiff(GetProms,names(grl))
fail
grl2 <- grl[intersect(names(grl),GetProms)]
#grl3 <- transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by="gene")[grl2]
promoter.seqs <- getPromoterSeq(grl2, Hsapiens, upstream=1000,
                                downstream=0)
#Spromoter.list <- sapply(grl,getPromoterSeq, Hsapiens, upstream=1000, dowstream=0)

### Search promoters for MYB binding sequence
# First convert 

#pfm.MYB.jaspar <- query(query (MotifDb, "mmusculus"),"MYB")[[3]]
#pcm.MYB.jaspar <- round(100 * pfm.MYB.jaspar)
#matchPWM(pcm.MYB.jaspar, unlist(promoter.seqs)[[1]], "90%")

print (class(promoter.seqs))
promoter.seqs <- unlist(promoter.seqs)
print (class(promoter.seqs))

pwm.hits <- sapply(promoter.seqs, 
                   function(pseq) 
                     matchPWM(pcm.MYB.jaspar, pseq, min.score="90%"))

Myb <- sapply(pwm.hits, length)
tbl.MYB <- data.frame(gene=genes, Myb)

### Sort sequences into MYB+ and MYB- piles

### Find over-represented motifs in MYB+ sequences ()

### Match over-represented motifs to other TF binding sites (MotIV??)




##   Views on a 1000-letter DNAString subject
## subject: TTGAGGAGTTGTCCACATACACATTAGTGTTG...AAAAAAAAGTGAAATACTGCGAAGAACAAAG
## views:
##     start end width
## [1]   620 626     7 [TGATAAG]
## [2]   637 643     7 [CGATAAG]