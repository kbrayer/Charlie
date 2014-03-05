#An attempt to look at TF binding sites in promoter regions of significant genes from RNA-Seq Data
#Charlie Brayer Feb 2014

#This script will take output from RNA seq data analysis, obtain promoter sequences,
#look for specific motifs and score promoters - somehow. Hopefully will eventually also
#scan promoters for unknown TF binding sites. 



library(MotifDb)
library(seqLogo)
library(motifStack)
library(Biostrings)
library(GenomicFeatures)
library(hgu133plus2.db)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

###Get promoter sequence for differentially expressed genes.
#Read in list of genes

SigGenes <- read.table("Significant_RNA_seq_genes.txt", header=TRUE, sep ='\t')

#Convert to Entrezgene Accession Numbers
Entrezgenes <- SigGenes[,3]
length(Entrezgenes)
head(Entrezgenes)
egREFSEQ <- toTable(org.Hs.egREFSEQ)
m <- match(SigGenes$ID, egREFSEQ$accession)
SigGenes$EntrezGene <- egREFSEQ$gene_id[m]
GetGenePromoters <- as.character(SigGenes [,8])

#Get Genome Coordinates for Genes of Interest
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
grl <- transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by="gene")
#Throw out genes with no genome location can't be found
fail <- setdiff(GetGenePromoters,names(grl))
fail
grl2 <- grl[intersect(names(grl),GetGenePromoters)]

#Get promoter sequences for genes of interest
#WARNING: this is slow...
promoter.seqs <- getPromoterSeq(grl2, Hsapiens, upstream=1000,
                                downstream=0)

### Search promoters for MYB binding sequence
#Load MYB binding sequence motif and convert from PFM to PCM

query(MotifDb, "MYB")
pfm.MYB.jaspar <- query (query (MotifDb, "mmusculus"),"MYB") [[3]]
seqLogo(pfm.MYB.jaspar)
pcm.MYB.jaspar <- round(100 * pfm.MYB.jaspar)

##Weird unlist stuff, no idea what it is or does...

print (class(promoter.seqs))
promoter.seqs <- unlist(promoter.seqs)
print (class(promoter.seqs))

#Scan promoters
###NOTE: Threshold may need to be changed, currently only finds 90% match


Myb.hits <- sapply(promoter.seqs, function (pseq) 
  matchPWM(pcm.MYB.jaspar, pseq, min.score = "80%", with.score=TRUE))

rMyb.hits <- sapply(promoter.seqs, function (pseq)
              matchPWM(reverseComplement(pcm.MYB.jaspar, pseq, min.score = "80%", with.score = TRUE)))



#Myb.hits <- sapply(promoter.seqs, function (pseq)
    #matchPWM(pcm.MYB.jaspar, pseq, min.score = "50%"))

myb.jaspar <- sapply(Myb.hits, length)
rmyb.jaspar <- sapply(rMyb.hits, length)

##convert Entrez accession back to Gene Symbols for making the output table readable

EntrezPromoters <- as.data.frame(names(promoter.seqs))
egSYMBOL <- toTable(org.Hs.egSYMBOL)
head(egSYMBOL)
m2 <- match(EntrezPromoters [,1], egSYMBOL$gene_id)
EntrezPromoters$symbol <- egSYMBOL$symbol[m2]

rEntrezPromoters <- as.data.frame(names(promoter.seqs))
egSYMBOL <- toTable(org.Hs.egSYMBOL)
head(egSYMBOL)
m2 <- match(rEntrezPromoters [,1], egSYMBOL$gene_id)
rEntrezPromoters$symbol <- egSYMBOL$symbol[m2]
rtbl.myb <- data.frame(gene=rEntrezPromoters$symbol, rmyb.jaspar) 

